# Changelog

## 0.2.0 - 2026-08-11

- Rename the PyPI distribution to `nanotime-ont`; keep the `nanotime` command.
- Add `nanotime inspect` for acquisition metadata, clock health, and yield checks.
- Add arbitrary elapsed-time checkpoints with `--checkpoints`.
- Add cumulative base-yield checkpoints with `--yield` and reached-time reporting.
- Add `nanotime timeline` for binned read and base throughput TSV output.
- Add zero-output split planning with `--dry-run` and BAM size estimates.
- Add the storage-efficient `nanotime manifest` disjoint workflow.
- Reject multiple `@RG DT` acquisition epochs unless `--allow-multiple-runs` is set.
- Warn about pre-origin reads, hour-long clock gaps, and multi-day timelines.
