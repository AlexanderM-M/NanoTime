# NanoTime

NanoTime reconstructs Oxford Nanopore sequencing timelines from BAM acquisition
metadata. It can inspect a run without writing output, create elapsed-time or
yield-based BAM checkpoints, and export a small throughput table for plotting.

```bash
nanotime inspect *.bam
nanotime split *.bam --interval 10m --until 90m --cumulative
```

NanoTime uses each read's `st:Z` timestamp and, by default, its `du:f` duration.
It does not assume that BAM chunk names, chunk counts, or file creation times
correspond to sequencing time.

## Installation

NanoTime requires Python 3.10 or newer. The PyPI distribution is named
`nanotime-ont` because the `nanotime` distribution belongs to an unrelated
project. The installed command remains `nanotime`.

```bash
python -m pip install nanotime-ont
nanotime --version
```

To install from source:

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

## Inspect before writing BAMs

`inspect` validates the clock, acquisition tags, run epochs, and expected yield
without creating any files:

```bash
nanotime inspect /data/run/pass/*.bam
```

It reports input and unique-read counts, the inferred UTC origin, first and last
read times, `st`/`du` coverage, supplementary records without `st`, yield at
common elapsed-time checkpoints, clock warnings, and a suggested split command.

## Split by elapsed time

Regular cumulative checkpoints:

```bash
nanotime split /data/run/pass/*.bam \
  --interval 10m \
  --until 90m \
  --cumulative \
  --prefix patient_01 \
  --output timed_bams/
```

This creates `patient_01_010min.bam` through `patient_01_090min.bam`. Each
output contains all reads that were complete before that checkpoint.

Researchers often need irregular checkpoints. Use a strictly increasing,
comma-separated list:

```bash
nanotime split *.bam \
  --checkpoints 5m,10m,15m,30m,60m,90m \
  --cumulative
```

This creates `sample_005min.bam`, `sample_010min.bam`,
`sample_015min.bam`, `sample_030min.bam`, `sample_060min.bam`, and
`sample_090min.bam`.

Disjoint output is the default for time checkpoints; `--disjoint` makes the
intent explicit:

```bash
nanotime split *.bam --interval 10m --until 30m --disjoint
```

The outputs are `sample_000-010min.bam`, `sample_010-020min.bam`, and
`sample_020-030min.bam`. If `--until` is omitted, NanoTime rounds up to the
first regular interval boundary after the final read.

Windows are half-open: `[start, end)`. A read available at exactly 10:00 is in
the 10–20 minute disjoint window and is not in the 10-minute cumulative BAM.

## Split by yield

Yield checkpoints compare experiments by the amount of sequence produced
rather than wall-clock time:

```bash
nanotime split *.bam --yield 100M,250M,500M,1G
```

Yield checkpoints are cumulative and create files such as
`sample_100Mb.bam` and `sample_1Gb.bam`. NanoTime includes the read that crosses
each threshold and prints when every checkpoint was reached. The same data is
recorded in `timeline_summary.tsv`.

## Export a throughput timeline

`timeline` aggregates unique reads and bases without creating BAMs:

```bash
nanotime timeline *.bam --bin 1m
nanotime timeline *.bam --bin 10m --output timeline.tsv
```

The TSV columns are:

```text
bin_start_seconds  bin_end_seconds  new_reads  cumulative_reads  new_bases  cumulative_bases
```

## Preview disk use

`--dry-run` performs timestamp scanning and assignment but writes no output
directory, BAM, index, summary, or manifest:

```bash
nanotime split *.bam \
  --interval 10m \
  --until 90m \
  --cumulative \
  --dry-run
```

It lists the estimated reads, alignment records, sequence yield, and compressed
BAM size for every planned output. The size estimate scales the input BAM bytes
by the number of assigned alignment records, so treat it as a planning estimate,
not an exact reservation.

## Storage-efficient checkpoint manifests

Cumulative BAMs repeatedly store early reads. The `manifest` workflow creates
disjoint BAMs plus `timeline_manifest.json` describing which files compose each
later checkpoint:

```bash
nanotime manifest *.bam --interval 10m --until 90m --output timed_bams/
```

For example, a 30-minute analysis consumes the 0–10, 10–20, and 20–30 minute
BAMs together. The sequence is stored once while the JSON manifest makes the
composition explicit for benchmarking pipelines.

## Clock origin and timestamp policy

Clock zero defaults to the earliest `@RG DT` header timestamp. If no read-group
start is available, NanoTime uses the earliest `st` value.

```bash
# Assign by pore entry time rather than read completion
nanotime split *.bam --interval 10m --timestamp start

# Make the earliest read clock zero
nanotime split *.bam --interval 10m --origin first-read

# Supply a known experiment start
nanotime split *.bam --interval 10m \
  --origin 2026-07-29T13:43:45.335073+02:00
```

Acquisition tags can be absent from secondary and supplementary alignments.
NanoTime maps timestamps by read ID first, then keeps every alignment for a read
in the same window.

## Multiple runs and suspicious clocks

NanoTime rejects distinct `@RG DT` acquisition epochs by default. This catches
the common mistake of passing BAMs from separate runs into one timeline. If the
merge is intentional, add `--allow-multiple-runs`; NanoTime still emits a
warning.

The scanner also warns when reads precede clock zero, a gap of at least one hour
exists between reads, or timestamps span at least one day. Run `inspect` first
when working with merged, resumed, or restarted acquisitions.

## Main split options

```text
--interval DURATION       Regular size: 30s, 10m, 1.5h, or 1d
--checkpoints D1,D2,...   Arbitrary elapsed-time checkpoints
--yield Y1,Y2,...         Cumulative base checkpoints: 100M, 1G, ...
--until DURATION          Stop regular intervals at elapsed time
--cumulative              Produce elapsed-time checkpoints
--disjoint                Produce non-overlapping windows (time default)
--dry-run                 Estimate outputs without writing files
--timestamp end|start     Use completion (default) or start time
--origin TIME             auto, first-read, or an ISO-8601 timestamp
--allow-multiple-runs     Permit distinct @RG DT acquisition epochs
--prefix NAME             Output prefix (default: sample)
--threads N               Sorting/indexing threads
--no-index                Do not create BAM indexes
--missing error|skip      Policy for reads with no timestamp
--fast                    Keep the timestamp map in RAM
--progress/--no-progress  Control progress output
--force                   Replace same-named NanoTime outputs
```

Run `nanotime split --help`, `nanotime inspect --help`, or
`nanotime timeline --help` for full command references. Input globs can be
expanded by the shell or quoted and expanded by NanoTime.

## Output guarantees and practical notes

- Output BAMs are coordinate-sorted and indexed by default.
- All alignments sharing a read ID stay together.
- Input BAMs must use compatible reference sequences and lengths.
- Existing same-named outputs are protected unless `--force` is supplied.
- Reads before clock zero and reads at or beyond an elapsed-time endpoint are
  excluded.
- Cumulative mode uses more time and disk than disjoint mode because early
  reads appear in every later BAM.
- The default temporary SQLite map keeps RAM use bounded; `--fast` trades RAM
  for speed.
- `timeline_summary.tsv` records UTC boundaries, unique reads, alignment
  records, and bases for each output.

## License

MIT
