# nanoTime

nanoTime splits Oxford Nanopore BAM files into exact experimental time
windows using each read's acquisition timestamp. It does not assume that BAM
chunk names, chunk counts, or file creation times correspond to sequencing
time.

```text
nanotime split *.bam \
  --interval 10m \
  --until 90m \
  --cumulative \
  --output timed_bams/
```

This creates:

```text
timed_bams/
├── sample_010min.bam
├── sample_010min.bam.bai
├── sample_020min.bam
├── sample_020min.bam.bai
├── ...
├── sample_090min.bam
├── sample_090min.bam.bai
└── timeline_summary.tsv
```

The 20-minute BAM contains all reads available during the first 20 minutes,
regardless of which MinKNOW/Dorado output chunk originally contained them.
This is useful for measuring when a real-time classifier, coverage estimate,
or variant call becomes reliable.

## How the clock works

Dorado/MinKNOW BAMs store a read's ISO-8601 start time in the `st:Z` tag and
its duration in seconds in `du:f`. nanoTime defaults to the read completion
time (`st + du`), because this is the earliest time the complete read could
have been available to downstream analysis.

Clock zero defaults to the earliest run start in the BAM `@RG DT` header. If
that metadata is absent, nanoTime uses the earliest read start. The policy
can be changed explicitly:

```text
# Assign by pore entry time rather than completion time
nanotime split *.bam --interval 10m --timestamp start

# Make the earliest read clock zero
nanotime split *.bam --interval 10m --origin first-read

# Supply a known experiment start
nanotime split *.bam --interval 10m \
  --origin 2026-07-29T13:43:45.335073+02:00
```

Windows are half-open: `[start, end)`. A read available at exactly 10:00
belongs to the 10–20 minute disjoint window and is not present in the
10-minute cumulative checkpoint.

Acquisition tags can be absent from secondary and supplementary alignment
records. nanoTime first maps timestamps by read ID, then keeps every
alignment for a read in the same time window.

## Installation

nanoTime requires Python 3.10 or newer. Installing it also installs
[pysam](https://pysam.readthedocs.io/), which bundles the required htslib
functionality.

Install from PyPI:

```bash
python -m pip install nanotime
```

Or install from GitHub:

```bash
git clone https://github.com/AlexanderM-M/NanoTime.git nanoTime
cd nanoTime
python -m pip install .
```

For development:

```bash
python -m pip install -e ".[test]"
pytest
```

## Usage

### Cumulative checkpoints

```bash
nanotime split /data/run/pass/*.bam \
  --interval 10m \
  --until 90m \
  --cumulative \
  --prefix patient_01 \
  --output timed_bams/
```

Each output contains all reads from clock zero through that checkpoint:

```text
patient_01_010min.bam
patient_01_020min.bam
...
patient_01_090min.bam
```

### Non-overlapping intervals

Disjoint output is the default; `--disjoint` makes the intent explicit:

```bash
nanotime split /data/run/pass/*.bam \
  --interval 10m \
  --until 90m \
  --disjoint \
  --prefix patient_01 \
  --output timed_bams/
```

This creates:

```text
patient_01_000-010min.bam
patient_01_010-020min.bam
patient_01_020-030min.bam
...
```

If `--until` is omitted, nanoTime rounds up to the first interval boundary
after the final read.

### Main options

```text
--interval DURATION     Required window size: 30s, 10m, 1.5h, or 1d
--until DURATION        Stop at elapsed time; defaults to the final read
--cumulative            Produce elapsed-time checkpoints
--disjoint              Produce non-overlapping windows (default)
--timestamp end|start   Use read completion (default) or read start
--origin TIME           auto, first-read, or an ISO-8601 timestamp
--prefix NAME           Output prefix (default: sample)
--threads N             Sorting/indexing threads
--no-index              Do not create BAM indexes
--missing error|skip    Policy for reads with no timestamp (default: error)
--fast                  Load read timestamps into RAM for faster splitting
--progress              Force progress output on non-interactive terminals
--no-progress           Hide progress output
--force                 Replace same-named nanoTime outputs
```

Run `nanotime split --help` for the complete command reference. Input globs
may be expanded by the shell or quoted and expanded by nanoTime.

## Timeline summary

`timeline_summary.tsv` records the exact UTC boundaries and yield in each
output:

```text
file	mode	start_seconds	end_seconds	start_time_utc	end_time_utc	reads	alignment_records	bases
sample_010min.bam	cumulative	0	600	2026-07-29T11:43:45.335073Z	2026-07-29T11:53:45.335073Z	...	...	...
```

`reads` and `bases` count unique read IDs. `alignment_records` includes
primary, secondary, and supplementary records written to the BAM.

## Output guarantees and practical notes

- Output BAMs are coordinate-sorted and indexed by default.
- All alignments sharing a read ID stay together.
- Input BAMs must use compatible reference sequences and lengths.
- Existing same-named outputs are protected unless `--force` is supplied.
- Reads before clock zero and reads at or after `--until` are excluded.
- Cumulative mode necessarily writes early reads into multiple files, so it
  uses more time and disk space than disjoint mode.
- By default, nanoTime uses a temporary on-disk SQLite map inside the output
  directory to keep memory usage bounded on large runs.
- `--fast` keeps the read timestamp map in RAM for faster splitting at the cost
  of higher memory usage.

## Why not split by BAM file?

MinKNOW and Dorado rotate output chunks for operational reasons such as batch
size and processing throughput. Eleven output BAMs therefore do not imply
eleven ten-minute periods. nanoTime treats the files as storage containers
and reconstructs experimental time from the read metadata itself.

## License

MIT
