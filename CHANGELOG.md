# Changelog

## 0.3.0 - 2026-08-12

- Add `nanotime plot` for a two-panel cumulative-yield and binned-throughput figure.
- Support PNG, SVG, and PDF plots with automatic elapsed-time and base-unit scaling.
- Add a documented example visualization and keep Matplotlib in the optional `plot` extra.
- Make PyPI publishing opt-in until a trusted publisher is configured.

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
