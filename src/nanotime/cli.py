"""Command-line interface for NanoTime."""

from __future__ import annotations

import argparse
import csv
import math
import shlex
import sys
from pathlib import Path
from typing import Sequence

from . import __version__
from .split import (
    MINUTE_NS,
    NanoTimeError,
    SplitConfig,
    dry_run_bams,
    execute_split,
    inspect_bams,
    timeline_bams,
)
from .timeutil import duration_ns, format_bytes, format_size, format_timestamp, parse_size


def _duration_argument(value: str) -> int:
    try:
        return duration_ns(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc


def _duration_list(value: str) -> tuple[int, ...]:
    try:
        items = tuple(duration_ns(item) for item in value.split(",") if item.strip())
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc
    if not items:
        raise argparse.ArgumentTypeError("provide at least one comma-separated checkpoint")
    if tuple(sorted(set(items))) != items:
        raise argparse.ArgumentTypeError("checkpoints must be unique and in increasing order")
    return items


def _yield_list(value: str) -> tuple[int, ...]:
    try:
        items = tuple(parse_size(item) for item in value.split(",") if item.strip())
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc
    if not items:
        raise argparse.ArgumentTypeError("provide at least one comma-separated yield checkpoint")
    if tuple(sorted(set(items))) != items:
        raise argparse.ArgumentTypeError("yield checkpoints must be unique and in increasing order")
    return items


def _add_inputs(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "inputs",
        nargs="+",
        metavar="BAM",
        help="input BAM file(s); quoted glob patterns are accepted",
    )


def _add_clock_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--timestamp",
        choices=("end", "start"),
        default="end",
        help="assign by read completion (st + du, default) or read start (st)",
    )
    parser.add_argument(
        "--origin",
        default="auto",
        metavar="TIME",
        help=(
            "clock zero: auto uses earliest @RG DT; alternatively use first-read "
            "or an ISO-8601 timestamp (default: auto)"
        ),
    )
    parser.add_argument(
        "--allow-multiple-runs",
        action="store_true",
        help="accept distinct acquisition epochs in the input (warning is still shown)",
    )
    parser.add_argument(
        "--fast",
        action="store_true",
        help="keep read timestamps in RAM instead of a temporary SQLite database",
    )


def _add_progress_options(parser: argparse.ArgumentParser) -> None:
    progress = parser.add_mutually_exclusive_group()
    progress.add_argument("--progress", action="store_true", dest="progress", help="show per-file progress updates")
    progress.add_argument("--no-progress", dest="progress", action="store_false", help="hide per-file progress updates")
    parser.set_defaults(progress=sys.stderr.isatty())


def _add_split_options(
    parser: argparse.ArgumentParser,
    *,
    include_mode: bool,
    include_yield: bool,
) -> None:
    _add_inputs(parser)
    checkpoints = parser.add_mutually_exclusive_group(required=True)
    checkpoints.add_argument(
        "--interval",
        type=_duration_argument,
        metavar="DURATION",
        help="regular window size, e.g. 10m, 30s, or 1.5h",
    )
    checkpoints.add_argument(
        "--checkpoints",
        type=_duration_list,
        metavar="D1,D2,...",
        help="arbitrary elapsed-time checkpoints, e.g. 5m,10m,30m,60m",
    )
    if include_yield:
        checkpoints.add_argument(
            "--yield",
            dest="yield_targets",
            type=_yield_list,
            metavar="Y1,Y2,...",
            help="cumulative base-yield checkpoints, e.g. 100M,500M,1G",
        )
    parser.add_argument(
        "--until",
        type=_duration_argument,
        metavar="DURATION",
        help="stop regular intervals at this elapsed time (default: after final read)",
    )
    if include_mode:
        mode = parser.add_mutually_exclusive_group()
        mode.add_argument(
            "--cumulative",
            action="store_true",
            help="each BAM contains all reads available from time zero to its endpoint",
        )
        mode.add_argument(
            "--disjoint",
            action="store_true",
            help="create non-overlapping intervals (the default for time checkpoints)",
        )
    parser.add_argument("--output", "-o", type=Path, default=Path("timed_bams"), help="output directory (default: timed_bams)")
    parser.add_argument("--prefix", default="sample", help="output filename prefix (default: sample)")
    _add_clock_options(parser)
    parser.add_argument("--threads", type=int, default=1, help="threads used to sort and index each BAM (default: 1)")
    parser.add_argument("--no-index", action="store_true", help="do not create .bam.bai index files")
    parser.add_argument(
        "--missing",
        choices=("error", "skip"),
        default="error",
        help="handling for reads with no st-tagged record (default: error)",
    )
    parser.add_argument("--force", action="store_true", help="replace same-named NanoTime outputs")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="scan and estimate every output without writing files",
    )
    _add_progress_options(parser)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nanotime",
        description="Inspect and reconstruct Oxford Nanopore acquisition timelines.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    split = subparsers.add_parser("split", help="create time- or yield-checkpoint BAMs")
    _add_split_options(split, include_mode=True, include_yield=True)

    manifest = subparsers.add_parser(
        "manifest",
        help="create storage-efficient disjoint BAMs and a composition manifest",
    )
    _add_split_options(manifest, include_mode=False, include_yield=False)

    inspect = subparsers.add_parser(
        "inspect",
        help="validate acquisition metadata and estimate elapsed-time yield",
    )
    _add_inputs(inspect)
    _add_clock_options(inspect)
    _add_progress_options(inspect)

    timeline = subparsers.add_parser(
        "timeline",
        help="write a small binned throughput table without creating BAMs",
    )
    _add_inputs(timeline)
    timeline.add_argument("--bin", required=True, type=_duration_argument, metavar="DURATION", help="timeline bin width, e.g. 1m or 10m")
    timeline.add_argument("--output", "-o", type=Path, help="write TSV to this path (default: standard output)")
    _add_clock_options(timeline)
    _add_progress_options(timeline)
    return parser


def _elapsed(value_ns: int | None, origin_ns: int) -> str:
    if value_ns is None:
        return "n/a"
    seconds = (value_ns - origin_ns) / 1_000_000_000
    sign = "+" if seconds >= 0 else "-"
    seconds = abs(seconds)
    hours, remainder = divmod(int(seconds), 3600)
    minutes, whole_seconds = divmod(remainder, 60)
    fraction = seconds - int(seconds)
    if fraction:
        return f"{sign}{hours:02d}:{minutes:02d}:{whole_seconds + fraction:06.3f}"
    return f"{sign}{hours:02d}:{minutes:02d}:{whole_seconds:02d}"


def _percent(numerator: int, denominator: int) -> str:
    return f"{numerator / denominator * 100:.2f}%" if denominator else "n/a"


def _print_warnings(warnings: Sequence[str], *, stream=sys.stderr) -> None:
    for warning in warnings:
        print(f"WARNING: {warning}", file=stream)


def _run_inspect(args: argparse.Namespace) -> int:
    report = inspect_bams(
        args.inputs,
        origin=args.origin,
        timestamp_mode=args.timestamp,
        allow_multiple_runs=args.allow_multiple_runs,
        fast=args.fast,
        progress=args.progress,
    )
    scan = report.scan
    print("NanoTime timeline inspection\n")
    print(f"{'Input BAMs':<24}{len(report.inputs):,}")
    print(f"{'Unique reads':<24}{scan.unique_reads:,}")
    print(f"{'Alignment records':<24}{scan.total_records:,}\n")
    print("Clock\n" + "─" * 42)
    print(f"{'Run origin':<24}{format_timestamp(scan.origin_ns)}")
    print(f"{'First read start':<24}{_elapsed(scan.min_start_ns, scan.origin_ns)}")
    print(f"{'First read complete':<24}{_elapsed(scan.min_complete_ns, scan.origin_ns)}")
    print(f"{'Last read complete':<24}{_elapsed(scan.max_complete_ns, scan.origin_ns)}\n")
    print("Acquisition metadata\n" + "─" * 42)
    print(f"{'reads with st':<24}{_percent(scan.reads_with_st, scan.unique_reads)}")
    print(f"{'reads with du':<24}{_percent(scan.reads_with_du, scan.unique_reads)}")
    print(
        f"{'supplementary w/o st':<24}"
        f"{_percent(scan.supplementary_without_st, scan.supplementary_records)}\n"
    )
    print("Yield\n" + "─" * 42)
    for minutes, reads, bases in report.yields:
        print(f"{minutes:>3} min{'':<16}{format_size(bases):>10}  ({reads:,} reads)")
    if scan.warnings:
        print("\nWarnings\n" + "─" * 42)
        _print_warnings(scan.warnings, stream=sys.stdout)
    duration_ns = max(1, scan.max_event_ns - scan.origin_ns)
    until_minutes = max(10, math.ceil(duration_ns / (10 * MINUTE_NS)) * 10)
    print("\nSuggested command:\n")
    print(
        f"nanotime split {shlex.join(str(path) for path in report.inputs)} "
        f"--interval 10m --until {until_minutes}m --cumulative"
    )
    return 0


def _write_timeline(args: argparse.Namespace) -> int:
    report = timeline_bams(
        args.inputs,
        args.bin,
        origin=args.origin,
        timestamp_mode=args.timestamp,
        allow_multiple_runs=args.allow_multiple_runs,
        fast=args.fast,
        progress=args.progress,
    )
    handle = sys.stdout
    close_handle = False
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        handle = args.output.open("w", newline="", encoding="utf-8")
        close_handle = True
    try:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["bin_start_seconds", "bin_end_seconds", "new_reads", "cumulative_reads", "new_bases", "cumulative_bases"]
        )
        for row in report.rows:
            writer.writerow(
                [
                    f"{row.start_ns / 1_000_000_000:g}",
                    f"{row.end_ns / 1_000_000_000:g}",
                    row.new_reads,
                    row.cumulative_reads,
                    row.new_bases,
                    row.cumulative_bases,
                ]
            )
    finally:
        if close_handle:
            handle.close()
    _print_warnings(report.scan.warnings)
    if args.output is not None:
        print(f"Wrote {len(report.rows)} timeline bin(s) to {args.output}")
    return 0


def _split_config(args: argparse.Namespace) -> SplitConfig:
    yield_targets = getattr(args, "yield_targets", None)
    cumulative = bool(getattr(args, "cumulative", False))
    if yield_targets and not getattr(args, "disjoint", False):
        cumulative = True
    manifest = args.command == "manifest"
    if manifest:
        cumulative = False
    return SplitConfig(
        inputs=args.inputs,
        output=args.output.resolve(),
        interval_ns=args.interval,
        until_ns=args.until,
        cumulative=cumulative,
        prefix=args.prefix,
        timestamp_mode=args.timestamp,
        origin=args.origin,
        threads=args.threads,
        index=not args.no_index,
        force=args.force,
        missing=args.missing,
        progress=args.progress,
        fast=args.fast,
        checkpoints_ns=args.checkpoints,
        yield_targets=yield_targets,
        allow_multiple_runs=args.allow_multiple_runs,
        manifest=manifest,
    )


def _run_split(args: argparse.Namespace) -> int:
    if getattr(args, "yield_targets", None) and getattr(args, "disjoint", False):
        raise NanoTimeError("yield checkpoints are cumulative; do not use --disjoint")
    config = _split_config(args)
    if args.dry_run:
        execution = dry_run_bams(config)
        print("Would create:\n")
        for window, reads, records, bases in zip(
            execution.windows,
            execution.read_counts,
            execution.record_counts,
            execution.base_counts,
        ):
            print(
                f"{window.path.name:<28}{reads:>12,} reads  "
                f"{records:>12,} records  ~{format_size(bases):>8}"
            )
        print(f"\nEstimated total BAM size: ~{format_bytes(execution.estimated_output_bytes)}")
        if execution.yield_targets:
            print("\ncheckpoint\treached_at")
            for target, reached in zip(execution.yield_targets, execution.reached_elapsed_ns):
                print(f"{format_size(target)}\t{_elapsed(execution.origin_ns + reached, execution.origin_ns)[1:]}")
        _print_warnings(execution.warnings)
        return 0

    print("Scanning timestamps, assigning alignments, and finalizing BAMs...", file=sys.stderr, flush=True)
    execution = execute_split(config)
    windows = execution.windows
    kind = "disjoint manifest" if config.manifest else (
        "cumulative" if config.cumulative else "disjoint"
    )
    extra = " and timeline_manifest.json" if config.manifest else ""
    print(
        f"Wrote {len(windows)} {kind} BAM file(s), timeline_summary.tsv{extra} to {config.output}"
    )
    if execution.yield_targets:
        print("\ncheckpoint\treached_at")
        for target, reached in zip(execution.yield_targets, execution.reached_elapsed_ns):
            print(f"{format_size(target)}\t{_elapsed(execution.origin_ns + reached, execution.origin_ns)[1:]}")
    _print_warnings(execution.warnings)
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "inspect":
            return _run_inspect(args)
        if args.command == "timeline":
            return _write_timeline(args)
        if args.command in {"split", "manifest"}:
            return _run_split(args)
    except (NanoTimeError, OSError, ValueError) as exc:
        print(f"nanotime: error: {exc}", file=sys.stderr)
        return 2
    return 0
