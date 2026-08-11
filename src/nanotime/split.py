"""BAM timeline analysis and splitting implementation."""

from __future__ import annotations

import bisect
import csv
import glob
import json
import math
import os
import sqlite3
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import pysam

from . import __version__
from .timeutil import format_size, format_timestamp, parse_timestamp, to_epoch_ns

SECOND_NS = 1_000_000_000
MINUTE_NS = 60 * SECOND_NS
HOUR_NS = 60 * MINUTE_NS
DAY_NS = 24 * HOUR_NS


class NanoTimeError(RuntimeError):
    """A user-facing NanoTime error."""


@dataclass(frozen=True)
class SplitConfig:
    inputs: Sequence[str]
    output: Path
    interval_ns: int | None
    until_ns: int | None
    cumulative: bool
    prefix: str = "sample"
    timestamp_mode: str = "end"
    origin: str = "auto"
    threads: int = 1
    index: bool = True
    force: bool = False
    missing: str = "error"
    progress: bool = True
    fast: bool = False
    checkpoints_ns: Sequence[int] | None = None
    yield_targets: Sequence[int] | None = None
    allow_multiple_runs: bool = False
    manifest: bool = False


@dataclass(frozen=True)
class RunEpoch:
    origin_ns: int
    files: tuple[Path, ...]


@dataclass(frozen=True)
class ClockGap:
    start_ns: int
    end_ns: int


@dataclass(frozen=True)
class ScanResult:
    origin_ns: int
    min_start_ns: int
    max_event_ns: int
    min_complete_ns: int | None
    max_complete_ns: int | None
    tagged_records: int
    total_records: int
    unique_reads: int
    reads_with_st: int
    reads_with_du: int
    supplementary_records: int
    supplementary_without_st: int
    header: dict
    read_metadata: dict[str, tuple[int, int, int, bool]] | None
    epochs: tuple[RunEpoch, ...]
    warnings: tuple[str, ...]


@dataclass(frozen=True)
class OutputWindow:
    index: int
    start_ns: int
    end_ns: int
    path: Path


@dataclass(frozen=True)
class SplitExecution:
    windows: tuple[OutputWindow, ...]
    origin_ns: int
    read_counts: tuple[int, ...]
    record_counts: tuple[int, ...]
    base_counts: tuple[int, ...]
    mode: str
    warnings: tuple[str, ...]
    estimated_output_bytes: int
    yield_targets: tuple[int, ...] = ()
    reached_elapsed_ns: tuple[int, ...] = ()


@dataclass(frozen=True)
class InspectionReport:
    inputs: tuple[Path, ...]
    scan: ScanResult
    yields: tuple[tuple[int, int, int], ...]


@dataclass(frozen=True)
class TimelineRow:
    index: int
    start_ns: int
    end_ns: int
    new_reads: int
    cumulative_reads: int
    new_bases: int
    cumulative_bases: int


@dataclass(frozen=True)
class TimelineReport:
    inputs: tuple[Path, ...]
    scan: ScanResult
    bin_ns: int
    rows: tuple[TimelineRow, ...]


def _progress_update(enabled: bool, phase: str, index: int, total: int, path: Path) -> None:
    if not enabled:
        return
    percent = (index / total * 100) if total else 100.0
    if sys.stderr.isatty():
        print(
            f"\r{phase:>9} {index}/{total} ({percent:>5.1f}%) {path.name:<48}",
            end="",
            file=sys.stderr,
            flush=True,
        )
    else:
        print(f"{phase}: {index}/{total}: {path}", file=sys.stderr)


def _progress_done(enabled: bool) -> None:
    if enabled and sys.stderr.isatty():
        print("", file=sys.stderr)


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


def _read_event_ns(
    record: pysam.AlignedSegment,
    mode: str,
    *,
    require_duration: bool,
) -> tuple[int, int, int, bool]:
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

    has_du = record.has_tag("du")
    duration = 0.0
    if has_du:
        try:
            duration = float(record.get_tag("du"))
        except (TypeError, ValueError) as exc:
            raise NanoTimeError(f"read {record.query_name!r} has invalid du tag") from exc
        if not math.isfinite(duration) or duration < 0:
            raise NanoTimeError(
                f"read {record.query_name!r} has invalid du tag {duration!r}"
            )
    elif mode == "end" and require_duration:
        raise NanoTimeError(
            f"read {record.query_name!r} has st but no du tag; "
            "use --timestamp start or supply BAMs containing du"
        )

    duration_ns = int(round(duration * SECOND_NS))
    event_ns = start_ns + duration_ns if mode == "end" else start_ns
    return start_ns, duration_ns, event_ns, has_du


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


def _cluster_epochs(origin_files: dict[int, set[Path]]) -> tuple[RunEpoch, ...]:
    """Cluster nearly-identical read-group timestamps into acquisition epochs."""
    epochs: list[RunEpoch] = []
    for origin_ns in sorted(origin_files):
        files = origin_files[origin_ns]
        if epochs and origin_ns - epochs[-1].origin_ns <= MINUTE_NS:
            previous = epochs[-1]
            epochs[-1] = RunEpoch(
                previous.origin_ns,
                tuple(sorted(set(previous.files).union(files))),
            )
        else:
            epochs.append(RunEpoch(origin_ns, tuple(sorted(files))))
    return tuple(epochs)


def _multiple_runs_message(epochs: Sequence[RunEpoch]) -> str:
    lines = ["Input contains multiple acquisition epochs."]
    for index, epoch in enumerate(epochs, start=1):
        lines.extend(
            [
                "",
                f"Run {index}:",
                f"  {format_timestamp(epoch.origin_ns)}",
                f"  {len(epoch.files)} BAM file(s)",
            ]
        )
    lines.extend(["", "Use --allow-multiple-runs if this is intentional."])
    return "\n".join(lines)


def _create_database(database: sqlite3.Connection) -> None:
    database.execute("PRAGMA synchronous = OFF")
    database.execute("PRAGMA journal_mode = MEMORY")
    database.execute("PRAGMA temp_store = FILE")
    database.execute(
        """
        CREATE TABLE reads (
            name TEXT PRIMARY KEY,
            start_ns INTEGER NOT NULL,
            duration_ns INTEGER NOT NULL,
            event_ns INTEGER NOT NULL,
            bases INTEGER NOT NULL,
            has_du INTEGER NOT NULL
        )
        """
    )
    database.execute("CREATE TABLE observed_names (name TEXT PRIMARY KEY)")


def _clock_warnings(
    *,
    database: sqlite3.Connection | None,
    fast_full: dict[str, tuple[int, int, int, int, bool]],
    origin_ns: int,
    epochs: Sequence[RunEpoch],
) -> tuple[str, ...]:
    if database is not None:
        before = database.execute(
            "SELECT COUNT(*), MIN(start_ns) FROM reads WHERE start_ns < ?", (origin_ns,)
        ).fetchone()
        event_values: Iterable[int] = (
            row[0] for row in database.execute("SELECT event_ns FROM reads ORDER BY event_ns")
        )
    else:
        starts = [item[0] for item in fast_full.values() if item[0] < origin_ns]
        before = (len(starts), min(starts) if starts else None)
        event_values = iter(sorted(item[2] for item in fast_full.values()))

    warnings: list[str] = []
    if before[0]:
        earliest_offset = (before[1] - origin_ns) / SECOND_NS
        warnings.append(
            f"{before[0]} read(s) have timestamps before clock zero; "
            f"earliest offset: {earliest_offset:.3f} s"
        )
    if len(epochs) > 1:
        warnings.append(f"input contains {len(epochs)} acquisition epochs")
    first_event: int | None = None
    previous_event: int | None = None
    last_event: int | None = None
    largest_gap: ClockGap | None = None
    for event_ns in event_values:
        if first_event is None:
            first_event = event_ns
        if previous_event is not None and event_ns - previous_event >= HOUR_NS and (
            largest_gap is None
            or event_ns - previous_event > largest_gap.end_ns - largest_gap.start_ns
        ):
            largest_gap = ClockGap(previous_event, event_ns)
        previous_event = event_ns
        last_event = event_ns
    if first_event is not None and last_event is not None:
        span_ns = last_event - first_event
        if span_ns >= DAY_NS:
            warnings.append(f"input timestamps span {span_ns / DAY_NS:.2f} days")
        if largest_gap is not None:
            warnings.append(
                f"{(largest_gap.end_ns - largest_gap.start_ns) / HOUR_NS:.2f}-hour gap "
                f"detected between reads ({format_timestamp(largest_gap.start_ns)} → "
                f"{format_timestamp(largest_gap.end_ns)})"
            )
    return tuple(warnings)


def scan_inputs(
    inputs: Sequence[Path],
    database: sqlite3.Connection | None,
    timestamp_mode: str,
    origin_option: str,
    *,
    use_fast: bool = False,
    show_progress: bool = False,
    allow_multiple_runs: bool = False,
    require_duration: bool = True,
) -> ScanResult:
    """Scan timestamps, acquisition epochs, clock health, and BAM headers."""
    if database is not None:
        _create_database(database)
    elif not use_fast:
        raise NanoTimeError("database is required unless --fast is enabled")

    first_header: dict | None = None
    signature: tuple | None = None
    read_groups: dict[str, dict] = {}
    header_origin_values: list[int] = []
    origin_files: dict[int, set[Path]] = {}
    min_start_ns: int | None = None
    max_event_ns: int | None = None
    max_complete_ns: int | None = None
    min_complete_ns: int | None = None
    tagged_records = 0
    total_records = 0
    supplementary_records = 0
    supplementary_without_st = 0

    fast_metadata: dict[str, tuple[int, int, int, bool]] | None = {} if use_fast else None
    fast_full: dict[str, tuple[int, int, int, int, bool]] = {} if use_fast else {}
    observed_names: set[str] = set()
    total_inputs = len(inputs)

    for file_index, path in enumerate(inputs, start=1):
        _progress_update(show_progress, "scan", file_index, total_inputs, path)
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
            current_origins = _header_origins(header_dict)
            header_origin_values.extend(current_origins)
            for value in current_origins:
                origin_files.setdefault(value, set()).add(path)
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
                for record in bam.fetch(until_eof=True):
                    total_records += 1
                    name = record.query_name
                    if use_fast:
                        observed_names.add(name)
                    else:
                        database.execute("INSERT OR IGNORE INTO observed_names VALUES (?)", (name,))
                    has_st = record.has_tag("st")
                    if record.is_supplementary:
                        supplementary_records += 1
                        if not has_st:
                            supplementary_without_st += 1
                    if not has_st:
                        continue
                    start_ns, duration_ns, event_ns, has_du = _read_event_ns(
                        record,
                        timestamp_mode,
                        require_duration=require_duration,
                    )
                    tagged_records += 1
                    bases = record.query_length or 0
                    values = (start_ns, duration_ns, event_ns, bases, has_du)
                    if use_fast and fast_metadata is not None:
                        existing = fast_full.get(name)
                        if existing is None:
                            fast_full[name] = values
                            fast_metadata[name] = (start_ns, event_ns, bases, has_du)
                        else:
                            if tuple(existing[:3]) != values[:3]:
                                raise NanoTimeError(f"read {name!r} has conflicting timestamp tags")
                            combined_du = existing[4] or has_du
                            if bases > existing[3] or combined_du != existing[4]:
                                kept_bases = max(bases, existing[3])
                                fast_full[name] = (*existing[:3], kept_bases, combined_du)
                                fast_metadata[name] = (existing[0], existing[2], kept_bases, combined_du)
                    else:
                        existing = database.execute(
                            "SELECT start_ns, duration_ns, event_ns, bases, has_du FROM reads WHERE name = ?",
                            (name,),
                        ).fetchone()
                        if existing is None:
                            database.execute("INSERT INTO reads VALUES (?, ?, ?, ?, ?, ?)", (name, *values))
                        elif tuple(existing[:3]) != values[:3]:
                            raise NanoTimeError(f"read {name!r} has conflicting timestamp tags")
                        elif bases > existing[3] or (has_du and not existing[4]):
                            database.execute(
                                "UPDATE reads SET bases = ?, has_du = ? WHERE name = ?",
                                (max(bases, existing[3]), int(bool(has_du or existing[4])), name),
                            )

                    min_start_ns = start_ns if min_start_ns is None else min(min_start_ns, start_ns)
                    max_event_ns = event_ns if max_event_ns is None else max(max_event_ns, event_ns)
                    if has_du:
                        complete_ns = start_ns + duration_ns
                        min_complete_ns = complete_ns if min_complete_ns is None else min(min_complete_ns, complete_ns)
                        max_complete_ns = complete_ns if max_complete_ns is None else max(max_complete_ns, complete_ns)
            except OSError as exc:
                raise NanoTimeError(f"failed while reading BAM {path}: {exc}") from exc
        if database is not None:
            database.commit()

    if min_start_ns is None or max_event_ns is None or first_header is None:
        raise NanoTimeError(
            "no reads with an st tag were found; NanoTime requires Dorado/MinKNOW acquisition timestamps"
        )

    epochs = _cluster_epochs(origin_files)
    if len(epochs) > 1 and not allow_multiple_runs:
        raise NanoTimeError(_multiple_runs_message(epochs))

    if origin_option == "auto":
        origin_ns = min(header_origin_values) if header_origin_values else min_start_ns
    elif origin_option == "first-read":
        origin_ns = min_start_ns
    else:
        try:
            origin_ns = to_epoch_ns(parse_timestamp(origin_option))
        except ValueError as exc:
            raise NanoTimeError(f"invalid --origin value: {exc}") from exc

    if use_fast and fast_metadata is not None:
        unique_reads = len(observed_names)
        reads_with_st = len(fast_metadata)
        reads_with_du = sum(1 for item in fast_metadata.values() if item[3])
    else:
        database.execute("CREATE INDEX reads_event_idx ON reads(event_ns)")
        database.commit()
        unique_reads = database.execute("SELECT COUNT(*) FROM observed_names").fetchone()[0]
        reads_with_st = database.execute("SELECT COUNT(*) FROM reads").fetchone()[0]
        reads_with_du = database.execute("SELECT COUNT(*) FROM reads WHERE has_du = 1").fetchone()[0]
    if read_groups:
        first_header["RG"] = list(read_groups.values())
    warnings = _clock_warnings(
        database=database,
        fast_full=fast_full,
        origin_ns=origin_ns,
        epochs=epochs,
    )
    return ScanResult(
        origin_ns=origin_ns,
        min_start_ns=min_start_ns,
        max_event_ns=max_event_ns,
        min_complete_ns=min_complete_ns,
        max_complete_ns=max_complete_ns,
        tagged_records=tagged_records,
        total_records=total_records,
        unique_reads=unique_reads,
        reads_with_st=reads_with_st,
        reads_with_du=reads_with_du,
        supplementary_records=supplementary_records,
        supplementary_without_st=supplementary_without_st,
        header=_add_program_record(first_header),
        read_metadata=fast_metadata,
        epochs=epochs,
        warnings=warnings,
    )


def _boundaries(interval_ns: int, until_ns: int) -> list[int]:
    boundaries = list(range(interval_ns, until_ns + 1, interval_ns))
    if not boundaries or boundaries[-1] != until_ns:
        boundaries.append(until_ns)
    return boundaries


def _duration_label(ns: int, width: int, unit: str) -> str:
    if unit == "min":
        return f"{ns // MINUTE_NS:0{width}d}"
    if ns % SECOND_NS == 0:
        return f"{ns // SECOND_NS:0{width}d}"
    text = f"{ns / SECOND_NS:.3f}".rstrip("0").rstrip(".")
    return text.replace(".", "p")


def make_windows(
    output: Path,
    prefix: str,
    boundaries: Sequence[int],
    cumulative: bool,
    labels: Sequence[str] | None = None,
) -> list[OutputWindow]:
    all_whole_minutes = all(value % MINUTE_NS == 0 for value in (0, *boundaries))
    maximum = max(boundaries) // (MINUTE_NS if all_whole_minutes else SECOND_NS)
    width = max(3, len(str(maximum)))
    unit = "min" if all_whole_minutes else "sec"
    windows: list[OutputWindow] = []
    previous = 0
    for index, end_ns in enumerate(boundaries):
        end_label = _duration_label(end_ns, width, unit)
        if labels is not None:
            filename = f"{prefix}_{labels[index]}.bam"
            start_ns = 0
        elif cumulative:
            filename = f"{prefix}_{end_label}{unit}.bam"
            start_ns = 0
        else:
            start_label = _duration_label(previous, width, unit)
            filename = f"{prefix}_{start_label}-{end_label}{unit}.bam"
            start_ns = previous
        windows.append(OutputWindow(index, start_ns, end_ns, output / filename))
        previous = end_ns
    return windows


def _check_outputs(windows: Sequence[OutputWindow], force: bool, manifest: bool) -> None:
    bam_targets = [window.path for window in windows]
    targets = bam_targets + [Path(str(path) + ".bai") for path in bam_targets]
    targets.append(windows[0].path.parent / "timeline_summary.tsv")
    if manifest:
        targets.append(windows[0].path.parent / "timeline_manifest.json")
    existing = [path for path in targets if path.exists()]
    if existing and not force:
        preview = ", ".join(str(path) for path in existing[:3])
        raise NanoTimeError(
            f"output already exists ({preview}); use --force to replace NanoTime outputs"
        )


def _lookup_event(database: sqlite3.Connection, name: str) -> tuple[int, int] | None:
    row = database.execute("SELECT event_ns, bases FROM reads WHERE name = ?", (name,)).fetchone()
    return (row[0], row[1]) if row else None


def _iter_metadata(
    database: sqlite3.Connection | None,
    read_metadata: dict[str, tuple[int, int, int, bool]] | None,
    *,
    ordered: bool = False,
) -> Iterable[tuple[int, int]]:
    if read_metadata is not None:
        values: Iterable[tuple[int, int]] = (
            (item[1], item[2]) for item in read_metadata.values()
        )
        return iter(sorted(values)) if ordered else values
    if database is None:
        raise NanoTimeError("missing read metadata")
    order = " ORDER BY event_ns" if ordered else ""
    return database.execute(f"SELECT event_ns, bases FROM reads{order}")


def _summarize_reads(
    metadata: Iterable[tuple[int, int]],
    origin_ns: int,
    boundaries: Sequence[int],
    cumulative: bool,
) -> tuple[list[int], list[int]]:
    counts = [0] * len(boundaries)
    bases = [0] * len(boundaries)
    for event_ns, read_bases in metadata:
        elapsed = event_ns - origin_ns
        index = bisect.bisect_right(boundaries, elapsed)
        if elapsed < 0 or index >= len(boundaries):
            continue
        counts[index] += 1
        bases[index] += read_bases
    if cumulative:
        for index in range(1, len(boundaries)):
            counts[index] += counts[index - 1]
            bases[index] += bases[index - 1]
    return counts, bases


def _yield_boundaries(
    metadata: Iterable[tuple[int, int]],
    origin_ns: int,
    targets: Sequence[int],
) -> tuple[list[int], list[int]]:
    boundaries: list[int] = []
    reached: list[int] = []
    cumulative = 0
    target_index = 0
    for event_ns, bases in metadata:
        if event_ns < origin_ns:
            continue
        cumulative += bases
        while target_index < len(targets) and cumulative >= targets[target_index]:
            elapsed = event_ns - origin_ns
            reached.append(elapsed)
            # Time windows are half-open. Advancing by 1 ns includes the read that
            # crossed the requested yield threshold.
            boundaries.append(elapsed + 1)
            target_index += 1
        if target_index == len(targets):
            break
    if target_index != len(targets):
        raise NanoTimeError(
            f"yield checkpoint {format_size(targets[target_index])} was not reached; "
            f"maximum observed yield is {format_size(cumulative)}"
        )
    return boundaries, reached


def _write_summary(
    path: Path,
    windows: Sequence[OutputWindow],
    origin_ns: int,
    read_counts: Sequence[int],
    record_counts: Sequence[int],
    base_counts: Sequence[int],
    mode: str,
    yield_targets: Sequence[int] = (),
    reached_elapsed_ns: Sequence[int] = (),
) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        header = [
            "file", "mode", "start_seconds", "end_seconds", "start_time_utc",
            "end_time_utc", "reads", "alignment_records", "bases",
        ]
        if yield_targets:
            header.extend(["yield_checkpoint_bases", "reached_at_seconds"])
        writer.writerow(header)
        for index, (window, reads, records, bases) in enumerate(
            zip(windows, read_counts, record_counts, base_counts)
        ):
            row: list[object] = [
                window.path.name,
                mode,
                f"{window.start_ns / SECOND_NS:g}",
                f"{window.end_ns / SECOND_NS:g}",
                format_timestamp(origin_ns + window.start_ns),
                format_timestamp(origin_ns + window.end_ns),
                reads,
                records,
                bases,
            ]
            if yield_targets:
                row.extend([yield_targets[index], f"{reached_elapsed_ns[index] / SECOND_NS:g}"])
            writer.writerow(row)


def _write_manifest(
    path: Path,
    windows: Sequence[OutputWindow],
    origin_ns: int,
    read_counts: Sequence[int],
    base_counts: Sequence[int],
) -> None:
    outputs = []
    accumulated: list[str] = []
    for window, reads, bases in zip(windows, read_counts, base_counts):
        accumulated.append(window.path.name)
        outputs.append(
            {
                "file": window.path.name,
                "start_seconds": window.start_ns / SECOND_NS,
                "end_seconds": window.end_ns / SECOND_NS,
                "start_time_utc": format_timestamp(origin_ns + window.start_ns),
                "end_time_utc": format_timestamp(origin_ns + window.end_ns),
                "reads": reads,
                "bases": bases,
                "combine_through_checkpoint": list(accumulated),
            }
        )
    payload = {
        "schema_version": 1,
        "mode": "disjoint",
        "origin_utc": format_timestamp(origin_ns),
        "description": "Combine listed disjoint BAMs through a checkpoint without duplicating earlier reads on disk.",
        "outputs": outputs,
    }
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def _resolve_windows(
    config: SplitConfig,
    scan: ScanResult,
    database: sqlite3.Connection | None,
) -> tuple[list[int], list[OutputWindow], list[int], tuple[int, ...]]:
    yield_targets = tuple(config.yield_targets or ())
    reached: list[int] = []
    if yield_targets:
        if not config.cumulative:
            raise NanoTimeError("yield checkpoints are cumulative; do not use --disjoint")
        boundaries, reached = _yield_boundaries(
            _iter_metadata(database, scan.read_metadata, ordered=True),
            scan.origin_ns,
            yield_targets,
        )
        labels = [f"{format_size(value, compact=True)}b" for value in yield_targets]
        windows = make_windows(config.output, config.prefix, boundaries, True, labels)
    elif config.checkpoints_ns:
        boundaries = list(config.checkpoints_ns)
        windows = make_windows(config.output, config.prefix, boundaries, config.cumulative)
    else:
        if config.interval_ns is None:
            raise NanoTimeError("one of --interval, --checkpoints, or --yield is required")
        if config.until_ns is None:
            elapsed = max(1, scan.max_event_ns - scan.origin_ns)
            until_ns = ((elapsed // config.interval_ns) + 1) * config.interval_ns
        else:
            until_ns = config.until_ns
        if until_ns <= 0:
            raise NanoTimeError("--until must be greater than zero")
        boundaries = _boundaries(config.interval_ns, until_ns)
        windows = make_windows(config.output, config.prefix, boundaries, config.cumulative)
    return boundaries, windows, reached, yield_targets


def _run_split(config: SplitConfig, *, write_outputs: bool) -> SplitExecution:
    inputs = expand_inputs(config.inputs)
    if config.threads < 1:
        raise NanoTimeError("--threads must be at least 1")
    if config.until_ns is not None and (config.checkpoints_ns or config.yield_targets):
        raise NanoTimeError("--until can only be used with --interval")
    if config.manifest and config.cumulative:
        raise NanoTimeError("a storage-efficient manifest requires disjoint output")
    if write_outputs:
        config.output.mkdir(parents=True, exist_ok=True)
    temp_parent = config.output if write_outputs else None

    with tempfile.TemporaryDirectory(prefix=".nanotime-", dir=temp_parent) as temp_name:
        temp_dir = Path(temp_name)
        database: sqlite3.Connection | None = None
        if not config.fast:
            database = sqlite3.connect(temp_dir / "reads.sqlite")
        try:
            scan_start = time.perf_counter()
            scan = scan_inputs(
                inputs,
                database,
                config.timestamp_mode,
                config.origin,
                use_fast=config.fast,
                show_progress=config.progress,
                allow_multiple_runs=config.allow_multiple_runs,
            )
            _progress_done(config.progress)
            if config.progress:
                print(
                    f"Scanned {scan.tagged_records} st-tagged records from {len(inputs)} BAM(s) "
                    f"in {time.perf_counter() - scan_start:.1f}s",
                    file=sys.stderr,
                )

            boundaries, windows, reached, yield_targets = _resolve_windows(config, scan, database)
            if write_outputs:
                _check_outputs(windows, config.force, config.manifest)

            raw_paths = [temp_dir / f"window-{item.index}.bam" for item in windows]
            writers = (
                [pysam.AlignmentFile(str(path), "wb", header=scan.header) for path in raw_paths]
                if write_outputs else []
            )
            record_counts = [0] * len(windows)
            missing_names: set[str] = set()
            assign_start = time.perf_counter()
            try:
                for input_index, input_path in enumerate(inputs, start=1):
                    _progress_update(config.progress, "assign", input_index, len(inputs), input_path)
                    with pysam.AlignmentFile(str(input_path), "rb", check_sq=False) as bam:
                        for record in bam.fetch(until_eof=True):
                            if scan.read_metadata is not None:
                                item = scan.read_metadata.get(record.query_name)
                                lookup = (item[1], item[2]) if item is not None else None
                            else:
                                lookup = _lookup_event(database, record.query_name)
                            if lookup is None:
                                if config.missing == "skip":
                                    missing_names.add(record.query_name)
                                    continue
                                raise NanoTimeError(
                                    f"read {record.query_name!r} has no st-tagged alignment in any input; "
                                    "use --missing skip to omit it"
                                )
                            event_ns, _ = lookup
                            elapsed = event_ns - scan.origin_ns
                            index = bisect.bisect_right(boundaries, elapsed)
                            if elapsed < 0 or index >= len(boundaries):
                                continue
                            targets: Iterable[int] = range(index, len(windows)) if config.cumulative else (index,)
                            for target in targets:
                                if write_outputs:
                                    writers[target].write(record)
                                record_counts[target] += 1
            finally:
                for writer in writers:
                    writer.close()
            _progress_done(config.progress)
            if config.progress:
                print(
                    f"Assigned reads to {len(windows)} window(s) in "
                    f"{time.perf_counter() - assign_start:.1f}s",
                    file=sys.stderr,
                )
                if missing_names:
                    print(f"Missing st tag for {len(missing_names)} query names", file=sys.stderr)

            read_counts, base_counts = _summarize_reads(
                _iter_metadata(database, scan.read_metadata),
                scan.origin_ns,
                boundaries,
                config.cumulative,
            )
            if write_outputs:
                final_start = time.perf_counter()
                for raw_path, window in zip(raw_paths, windows):
                    try:
                        pysam.sort("-@", str(config.threads), "-o", str(window.path), str(raw_path))
                        if config.index:
                            pysam.index("-@", str(config.threads), str(window.path))
                    except pysam.SamtoolsError as exc:
                        raise NanoTimeError(
                            f"samtools failed while finalizing {window.path.name}: {exc}"
                        ) from exc
                if config.progress:
                    print(
                        f"Finalized {len(windows)} BAM(s) in {time.perf_counter() - final_start:.1f}s",
                        file=sys.stderr,
                    )
                mode = "yield-cumulative" if yield_targets else (
                    "cumulative" if config.cumulative else "disjoint"
                )
                _write_summary(
                    config.output / "timeline_summary.tsv",
                    windows,
                    scan.origin_ns,
                    read_counts,
                    record_counts,
                    base_counts,
                    mode,
                    yield_targets,
                    reached,
                )
                if config.manifest:
                    _write_manifest(
                        config.output / "timeline_manifest.json",
                        windows,
                        scan.origin_ns,
                        read_counts,
                        base_counts,
                    )
            mode = "yield-cumulative" if yield_targets else (
                "cumulative" if config.cumulative else "disjoint"
            )
            input_bytes = sum(path.stat().st_size for path in inputs)
            estimated = round(
                input_bytes * sum(record_counts) / scan.total_records
            ) if scan.total_records else 0
            return SplitExecution(
                tuple(windows),
                scan.origin_ns,
                tuple(read_counts),
                tuple(record_counts),
                tuple(base_counts),
                mode,
                scan.warnings,
                estimated,
                yield_targets,
                tuple(reached),
            )
        finally:
            if database is not None:
                database.close()


def split_bams(config: SplitConfig) -> list[OutputWindow]:
    """Split BAM inputs according to acquisition timestamps."""
    return list(execute_split(config).windows)


def execute_split(config: SplitConfig) -> SplitExecution:
    """Split BAM inputs and return the completed output statistics."""
    return _run_split(config, write_outputs=True)


def dry_run_bams(config: SplitConfig) -> SplitExecution:
    """Calculate split outputs and sizes without creating any output files."""
    return _run_split(config, write_outputs=False)


def inspect_bams(
    patterns: Sequence[str],
    *,
    origin: str = "auto",
    timestamp_mode: str = "end",
    allow_multiple_runs: bool = False,
    fast: bool = False,
    progress: bool = False,
) -> InspectionReport:
    """Inspect acquisition metadata and common elapsed-time yields."""
    inputs = expand_inputs(patterns)
    with tempfile.TemporaryDirectory(prefix="nanotime-inspect-") as temp_name:
        database = None if fast else sqlite3.connect(Path(temp_name) / "reads.sqlite")
        try:
            scan = scan_inputs(
                inputs,
                database,
                timestamp_mode,
                origin,
                use_fast=fast,
                show_progress=progress,
                allow_multiple_runs=allow_multiple_runs,
                require_duration=False,
            )
            yields: list[tuple[int, int, int]] = []
            for checkpoint in (10, 20, 30, 60, 90):
                boundary = checkpoint * MINUTE_NS
                if scan.read_metadata is not None:
                    chosen = [
                        (item[1], item[2])
                        for item in scan.read_metadata.values()
                        if 0 <= item[1] - scan.origin_ns < boundary
                    ]
                    read_count = len(chosen)
                    base_count = sum(item[1] for item in chosen)
                else:
                    read_count, base_count = database.execute(
                        "SELECT COUNT(*), COALESCE(SUM(bases), 0) FROM reads "
                        "WHERE event_ns >= ? AND event_ns < ?",
                        (scan.origin_ns, scan.origin_ns + boundary),
                    ).fetchone()
                yields.append((checkpoint, read_count, base_count))
            return InspectionReport(tuple(inputs), scan, tuple(yields))
        finally:
            if database is not None:
                database.close()


def timeline_bams(
    patterns: Sequence[str],
    bin_ns: int,
    *,
    origin: str = "auto",
    timestamp_mode: str = "end",
    allow_multiple_runs: bool = False,
    fast: bool = False,
    progress: bool = False,
) -> TimelineReport:
    """Aggregate unique read completions and yield into elapsed-time bins."""
    inputs = expand_inputs(patterns)
    with tempfile.TemporaryDirectory(prefix="nanotime-timeline-") as temp_name:
        database = None if fast else sqlite3.connect(Path(temp_name) / "reads.sqlite")
        try:
            scan = scan_inputs(
                inputs,
                database,
                timestamp_mode,
                origin,
                use_fast=fast,
                show_progress=progress,
                allow_multiple_runs=allow_multiple_runs,
            )
            elapsed_max = max(0, scan.max_event_ns - scan.origin_ns)
            bin_count = max(1, elapsed_max // bin_ns + 1)
            if bin_count > 1_000_000:
                raise NanoTimeError(
                    f"--bin would create {bin_count:,} rows; choose a bin width of at least "
                    f"{math.ceil(elapsed_max / 1_000_000 / SECOND_NS):g}s"
                )
            counts = [0] * bin_count
            bases = [0] * bin_count
            for event_ns, read_bases in _iter_metadata(database, scan.read_metadata):
                elapsed = event_ns - scan.origin_ns
                if elapsed < 0:
                    continue
                index = min(elapsed // bin_ns, bin_count - 1)
                counts[index] += 1
                bases[index] += read_bases
            rows: list[TimelineRow] = []
            cumulative_reads = 0
            cumulative_bases = 0
            for index, (new_reads, new_bases) in enumerate(zip(counts, bases)):
                cumulative_reads += new_reads
                cumulative_bases += new_bases
                rows.append(
                    TimelineRow(
                        index,
                        index * bin_ns,
                        (index + 1) * bin_ns,
                        new_reads,
                        cumulative_reads,
                        new_bases,
                        cumulative_bases,
                    )
                )
            return TimelineReport(tuple(inputs), scan, bin_ns, tuple(rows))
        finally:
            if database is not None:
                database.close()
