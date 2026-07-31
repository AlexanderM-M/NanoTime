"""Core BAM splitting implementation."""

from __future__ import annotations

import bisect
import csv
import glob
import math
import os
import sqlite3
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import pysam

from . import __version__
from .timeutil import format_timestamp, parse_timestamp, to_epoch_ns


class NanoTimeError(RuntimeError):
    """A user-facing nanoTime error."""


@dataclass(frozen=True)
class SplitConfig:
    inputs: Sequence[str]
    output: Path
    interval_ns: int
    until_ns: int | None
    cumulative: bool
    prefix: str = "sample"
    timestamp_mode: str = "end"
    origin: str = "auto"
    threads: int = 1
    index: bool = True
    force: bool = False
    missing: str = "error"


@dataclass(frozen=True)
class ScanResult:
    origin_ns: int
    min_start_ns: int
    max_event_ns: int
    tagged_records: int
    unique_reads: int
    header: dict


@dataclass(frozen=True)
class OutputWindow:
    index: int
    start_ns: int
    end_ns: int
    path: Path


def expand_inputs(patterns: Sequence[str]) -> list[Path]:
    """Expand shell-style input patterns and return unique BAM paths."""
    paths: list[Path] = []
    seen: set[Path] = set()
    for pattern in patterns:
        matches = glob.glob(os.path.expanduser(pattern))
        if not matches and Path(pattern).exists():
            matches = [pattern]
        if not matches:
            raise NanoTimeError(f"input did not match any files: {pattern}")
        for match in matches:
            path = Path(match).resolve()
            if path.suffix.lower() != ".bam":
                raise NanoTimeError(f"input is not a BAM file: {path}")
            if path not in seen:
                paths.append(path)
                seen.add(path)
    if not paths:
        raise NanoTimeError("no input BAM files were provided")
    return paths


def _header_signature(header: pysam.AlignmentHeader) -> tuple:
    return tuple(zip(header.references, header.lengths))


def _header_origins(header_dict: dict) -> list[int]:
    origins: list[int] = []
    for read_group in header_dict.get("RG", []):
        value = read_group.get("DT")
        if value:
            try:
                origins.append(to_epoch_ns(parse_timestamp(value)))
            except ValueError as exc:
                raise NanoTimeError(f"invalid @RG DT value {value!r}: {exc}") from exc
    return origins


def _read_event_ns(record: pysam.AlignedSegment, mode: str) -> tuple[int, int, int]:
    try:
        start_value = record.get_tag("st")
    except KeyError as exc:
        raise NanoTimeError("record has no st tag") from exc
    try:
        start_ns = to_epoch_ns(parse_timestamp(str(start_value)))
    except ValueError as exc:
        raise NanoTimeError(
            f"read {record.query_name!r} has invalid st tag {start_value!r}: {exc}"
        ) from exc

    duration = 0.0
    if record.has_tag("du"):
        try:
            duration = float(record.get_tag("du"))
        except (TypeError, ValueError) as exc:
            raise NanoTimeError(
                f"read {record.query_name!r} has invalid du tag"
            ) from exc
        if not math.isfinite(duration) or duration < 0:
            raise NanoTimeError(
                f"read {record.query_name!r} has invalid du tag {duration!r}"
            )
    elif mode == "end":
        raise NanoTimeError(
            f"read {record.query_name!r} has st but no du tag; "
            "use --timestamp start or supply BAMs containing du"
        )

    duration_ns = int(round(duration * 1_000_000_000))
    event_ns = start_ns + duration_ns if mode == "end" else start_ns
    return start_ns, duration_ns, event_ns


def _add_program_record(header: dict) -> dict:
    result = dict(header)
    programs = [dict(item) for item in result.get("PG", [])]
    ids = {item.get("ID") for item in programs}
    program_id = "nanotime"
    counter = 1
    while program_id in ids:
        counter += 1
        program_id = f"nanotime.{counter}"
    programs.append(
        {
            "ID": program_id,
            "PN": "nanotime",
            "VN": __version__,
            "DS": "split by Oxford Nanopore acquisition time",
        }
    )
    result["PG"] = programs
    return result


def scan_inputs(
    inputs: Sequence[Path],
    database: sqlite3.Connection,
    timestamp_mode: str,
    origin_option: str,
) -> ScanResult:
    """Scan timestamps into SQLite and validate the BAM headers."""
    database.execute(
        """
        CREATE TABLE reads (
            name TEXT PRIMARY KEY,
            start_ns INTEGER NOT NULL,
            duration_ns INTEGER NOT NULL,
            event_ns INTEGER NOT NULL,
            bases INTEGER NOT NULL
        )
        """
    )
    first_header: dict | None = None
    signature: tuple | None = None
    read_groups: dict[str, dict] = {}
    header_origin_values: list[int] = []
    min_start_ns: int | None = None
    max_event_ns: int | None = None
    tagged_records = 0

    for path in inputs:
        try:
            bam = pysam.AlignmentFile(str(path), "rb", check_sq=False)
        except (OSError, ValueError) as exc:
            raise NanoTimeError(f"cannot open BAM {path}: {exc}") from exc
        with bam:
            header_dict = bam.header.to_dict()
            current_signature = _header_signature(bam.header)
            if signature is None:
                signature = current_signature
                first_header = header_dict
            elif current_signature != signature:
                raise NanoTimeError(
                    f"BAM reference headers are incompatible: {inputs[0]} and {path}"
                )
            header_origin_values.extend(_header_origins(header_dict))
            for read_group in header_dict.get("RG", []):
                read_group_id = read_group.get("ID")
                if not read_group_id:
                    raise NanoTimeError(f"BAM has an @RG record without ID: {path}")
                existing_group = read_groups.get(read_group_id)
                if existing_group is None:
                    read_groups[read_group_id] = dict(read_group)
                elif existing_group != read_group:
                    raise NanoTimeError(
                        f"read group {read_group_id!r} differs between input BAMs"
                    )

            try:
                records = bam.fetch(until_eof=True)
                for record in records:
                    if not record.has_tag("st"):
                        continue
                    start_ns, read_duration_ns, event_ns = _read_event_ns(
                        record, timestamp_mode
                    )
                    tagged_records += 1
                    bases = record.query_length or 0
                    existing = database.execute(
                        "SELECT start_ns, duration_ns, event_ns, bases "
                        "FROM reads WHERE name = ?",
                        (record.query_name,),
                    ).fetchone()
                    values = (start_ns, read_duration_ns, event_ns, bases)
                    if existing is None:
                        database.execute(
                            "INSERT INTO reads VALUES (?, ?, ?, ?, ?)",
                            (record.query_name, *values),
                        )
                    elif tuple(existing[:3]) != values[:3]:
                        raise NanoTimeError(
                            f"read {record.query_name!r} has conflicting timestamp tags"
                        )
                    elif bases > existing[3]:
                        database.execute(
                            "UPDATE reads SET bases = ? WHERE name = ?",
                            (bases, record.query_name),
                        )
                    min_start_ns = (
                        start_ns if min_start_ns is None else min(min_start_ns, start_ns)
                    )
                    max_event_ns = (
                        event_ns if max_event_ns is None else max(max_event_ns, event_ns)
                    )
            except OSError as exc:
                raise NanoTimeError(f"failed while reading BAM {path}: {exc}") from exc
        database.commit()

    if min_start_ns is None or max_event_ns is None or first_header is None:
        raise NanoTimeError(
            "no reads with an st tag were found; nanoTime requires Dorado/MinKNOW "
            "acquisition timestamps"
        )

    if origin_option == "auto":
        origin_ns = min(header_origin_values) if header_origin_values else min_start_ns
    elif origin_option == "first-read":
        origin_ns = min_start_ns
    else:
        try:
            origin_ns = to_epoch_ns(parse_timestamp(origin_option))
        except ValueError as exc:
            raise NanoTimeError(f"invalid --origin value: {exc}") from exc

    unique_reads = database.execute("SELECT COUNT(*) FROM reads").fetchone()[0]
    if read_groups:
        first_header["RG"] = list(read_groups.values())
    return ScanResult(
        origin_ns=origin_ns,
        min_start_ns=min_start_ns,
        max_event_ns=max_event_ns,
        tagged_records=tagged_records,
        unique_reads=unique_reads,
        header=_add_program_record(first_header),
    )


def _boundaries(interval_ns: int, until_ns: int) -> list[int]:
    boundaries = list(range(interval_ns, until_ns + 1, interval_ns))
    if not boundaries or boundaries[-1] != until_ns:
        boundaries.append(until_ns)
    return boundaries


def _duration_label(ns: int, width: int, unit: str) -> str:
    if unit == "min":
        return f"{ns // 60_000_000_000:0{width}d}"
    if ns % 1_000_000_000 == 0:
        return f"{ns // 1_000_000_000:0{width}d}"
    text = f"{ns / 1_000_000_000:.3f}".rstrip("0").rstrip(".")
    return text.replace(".", "p")


def make_windows(
    output: Path,
    prefix: str,
    boundaries: Sequence[int],
    cumulative: bool,
) -> list[OutputWindow]:
    all_whole_minutes = all(
        value % 60_000_000_000 == 0 for value in (0, *boundaries)
    )
    maximum = (
        max(boundaries) // 60_000_000_000
        if all_whole_minutes
        else max(boundaries) // 1_000_000_000
    )
    width = max(3, len(str(maximum)))
    unit = "min" if all_whole_minutes else "sec"
    windows: list[OutputWindow] = []
    previous = 0
    for index, end_ns in enumerate(boundaries):
        end_label = _duration_label(end_ns, width, unit)
        if cumulative:
            filename = f"{prefix}_{end_label}{unit}.bam"
            start_ns = 0
        else:
            start_label = _duration_label(previous, width, unit)
            filename = f"{prefix}_{start_label}-{end_label}{unit}.bam"
            start_ns = previous
        windows.append(OutputWindow(index, start_ns, end_ns, output / filename))
        previous = end_ns
    return windows


def _check_outputs(windows: Sequence[OutputWindow], force: bool) -> None:
    bam_targets = [window.path for window in windows]
    targets = bam_targets + [Path(str(path) + ".bai") for path in bam_targets]
    targets.append(windows[0].path.parent / "timeline_summary.tsv")
    existing = [path for path in targets if path.exists()]
    if existing and not force:
        preview = ", ".join(str(path) for path in existing[:3])
        raise NanoTimeError(
            f"output already exists ({preview}); use --force to replace nanoTime outputs"
        )


def _lookup_event(
    database: sqlite3.Connection, name: str
) -> tuple[int, int] | None:
    row = database.execute(
        "SELECT event_ns, bases FROM reads WHERE name = ?", (name,)
    ).fetchone()
    return (row[0], row[1]) if row else None


def _write_summary(
    path: Path,
    windows: Sequence[OutputWindow],
    origin_ns: int,
    read_counts: Sequence[int],
    record_counts: Sequence[int],
    base_counts: Sequence[int],
    mode: str,
) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "file",
                "mode",
                "start_seconds",
                "end_seconds",
                "start_time_utc",
                "end_time_utc",
                "reads",
                "alignment_records",
                "bases",
            ]
        )
        for window, reads, records, bases in zip(
            windows, read_counts, record_counts, base_counts
        ):
            writer.writerow(
                [
                    window.path.name,
                    mode,
                    f"{window.start_ns / 1_000_000_000:g}",
                    f"{window.end_ns / 1_000_000_000:g}",
                    format_timestamp(origin_ns + window.start_ns),
                    format_timestamp(origin_ns + window.end_ns),
                    reads,
                    records,
                    bases,
                ]
            )


def _summarize_reads(
    database: sqlite3.Connection,
    origin_ns: int,
    boundaries: Sequence[int],
    cumulative: bool,
) -> tuple[list[int], list[int]]:
    counts = [0] * len(boundaries)
    bases = [0] * len(boundaries)
    for event_ns, read_bases in database.execute("SELECT event_ns, bases FROM reads"):
        elapsed = event_ns - origin_ns
        index = bisect.bisect_right(boundaries, elapsed)
        if elapsed < 0 or index >= len(boundaries):
            continue
        if cumulative:
            for target in range(index, len(boundaries)):
                counts[target] += 1
                bases[target] += read_bases
        else:
            counts[index] += 1
            bases[index] += read_bases
    return counts, bases


def split_bams(config: SplitConfig) -> list[OutputWindow]:
    """Split BAM inputs according to acquisition timestamps."""
    inputs = expand_inputs(config.inputs)
    config.output.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix=".nanotime-", dir=config.output) as temp_name:
        temp_dir = Path(temp_name)
        database = sqlite3.connect(temp_dir / "reads.sqlite")
        try:
            scan = scan_inputs(
                inputs, database, config.timestamp_mode, config.origin
            )
            if config.until_ns is None:
                elapsed = max(1, scan.max_event_ns - scan.origin_ns)
                # Windows are half-open, so an event exactly on a boundary belongs
                # in the next window.
                until_ns = (
                    (elapsed // config.interval_ns) + 1
                ) * config.interval_ns
            else:
                until_ns = config.until_ns
            if until_ns <= 0:
                raise NanoTimeError("--until must be greater than zero")

            boundaries = _boundaries(config.interval_ns, until_ns)
            windows = make_windows(
                config.output, config.prefix, boundaries, config.cumulative
            )
            _check_outputs(windows, config.force)

            raw_paths = [temp_dir / f"window-{item.index}.bam" for item in windows]
            writers = [
                pysam.AlignmentFile(str(path), "wb", header=scan.header)
                for path in raw_paths
            ]
            record_counts = [0] * len(windows)
            missing_names: set[str] = set()
            try:
                for input_path in inputs:
                    with pysam.AlignmentFile(
                        str(input_path), "rb", check_sq=False
                    ) as bam:
                        for record in bam.fetch(until_eof=True):
                            lookup = _lookup_event(database, record.query_name)
                            if lookup is None:
                                if config.missing == "skip":
                                    missing_names.add(record.query_name)
                                    continue
                                raise NanoTimeError(
                                    f"read {record.query_name!r} has no st-tagged "
                                    "alignment in any input; use --missing skip to omit it"
                                )
                            event_ns, _ = lookup
                            elapsed = event_ns - scan.origin_ns
                            index = bisect.bisect_right(boundaries, elapsed)
                            if elapsed < 0 or index >= len(boundaries):
                                continue
                            targets: Iterable[int]
                            if config.cumulative:
                                targets = range(index, len(windows))
                            else:
                                targets = (index,)
                            for target in targets:
                                writers[target].write(record)
                                record_counts[target] += 1
            finally:
                for writer in writers:
                    writer.close()

            for raw_path, window in zip(raw_paths, windows):
                try:
                    pysam.sort(
                        "-@",
                        str(config.threads),
                        "-o",
                        str(window.path),
                        str(raw_path),
                    )
                    if config.index:
                        pysam.index("-@", str(config.threads), str(window.path))
                except pysam.SamtoolsError as exc:
                    raise NanoTimeError(
                        f"samtools failed while finalizing {window.path.name}: {exc}"
                    ) from exc

            read_counts, base_counts = _summarize_reads(
                database, scan.origin_ns, boundaries, config.cumulative
            )
            _write_summary(
                config.output / "timeline_summary.tsv",
                windows,
                scan.origin_ns,
                read_counts,
                record_counts,
                base_counts,
                "cumulative" if config.cumulative else "disjoint",
            )
            return windows
        finally:
            database.close()
