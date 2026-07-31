"""Command-line interface for nanoTime."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .split import NanoTimeError, SplitConfig, split_bams
from .timeutil import duration_ns


def _duration_argument(value: str) -> int:
    try:
        return duration_ns(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nanotime",
        description=(
            "Split Oxford Nanopore BAM files by their true acquisition timestamps."
        ),
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    split = subparsers.add_parser(
        "split",
        help="create cumulative or non-overlapping time-window BAMs",
        description=(
            "Assign every read using its Dorado st tag and, by default, its du "
            "duration. All alignments of a read are kept together."
        ),
    )
    split.add_argument(
        "inputs",
        nargs="+",
        metavar="BAM",
        help="input BAM file(s); quoted glob patterns are accepted",
    )
    split.add_argument(
        "--interval",
        required=True,
        type=_duration_argument,
        metavar="DURATION",
        help="window size, e.g. 10m, 30s, or 1.5h",
    )
    split.add_argument(
        "--until",
        type=_duration_argument,
        metavar="DURATION",
        help="stop at this elapsed time (default: round up after the last read)",
    )
    mode = split.add_mutually_exclusive_group()
    mode.add_argument(
        "--cumulative",
        action="store_true",
        help="each BAM contains all reads available from time zero to its endpoint",
    )
    mode.add_argument(
        "--disjoint",
        action="store_true",
        help="create non-overlapping intervals (the default)",
    )
    split.add_argument(
        "--output",
        "-o",
        type=Path,
        default=Path("timed_bams"),
        help="output directory (default: timed_bams)",
    )
    split.add_argument(
        "--prefix",
        default="sample",
        help="output filename prefix (default: sample)",
    )
    split.add_argument(
        "--timestamp",
        choices=("end", "start"),
        default="end",
        help=(
            "assign by read completion (st + du, default) or read start (st)"
        ),
    )
    split.add_argument(
        "--origin",
        default="auto",
        metavar="TIME",
        help=(
            "clock zero: auto uses earliest @RG DT; alternatively use first-read "
            "or an ISO-8601 timestamp (default: auto)"
        ),
    )
    split.add_argument(
        "--threads",
        type=int,
        default=1,
        help="threads used to sort and index each BAM (default: 1)",
    )
    split.add_argument(
        "--no-index",
        action="store_true",
        help="do not create .bam.bai index files",
    )
    split.add_argument(
        "--missing",
        choices=("error", "skip"),
        default="error",
        help="handling for reads with no st-tagged record (default: error)",
    )
    split.add_argument(
        "--force",
        action="store_true",
        help="replace nanoTime output files with the same names",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "split":
        if args.threads < 1:
            parser.error("--threads must be at least 1")
        config = SplitConfig(
            inputs=args.inputs,
            output=args.output.resolve(),
            interval_ns=args.interval,
            until_ns=args.until,
            cumulative=args.cumulative,
            prefix=args.prefix,
            timestamp_mode=args.timestamp,
            origin=args.origin,
            threads=args.threads,
            index=not args.no_index,
            force=args.force,
            missing=args.missing,
        )
        print(
            "Scanning timestamps, assigning alignments, and finalizing BAMs...",
            file=sys.stderr,
            flush=True,
        )
        try:
            windows = split_bams(config)
        except (NanoTimeError, OSError, ValueError) as exc:
            print(f"nanotime: error: {exc}", file=sys.stderr)
            return 2
        print(
            f"Wrote {len(windows)} "
            f"{'cumulative' if config.cumulative else 'disjoint'} BAM file(s) "
            f"and timeline_summary.tsv to {config.output}"
        )
        return 0
    return 0
