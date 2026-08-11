from pathlib import Path

import pysam
import pytest

import json

from nanotime.split import (
    NanoTimeError,
    SplitConfig,
    dry_run_bams,
    inspect_bams,
    split_bams,
    timeline_bams,
)


def _record(
    name: str,
    position: int,
    start: str | None,
    duration: float | None,
    flag: int = 0,
) -> pysam.AlignedSegment:
    record = pysam.AlignedSegment()
    record.query_name = name
    record.query_sequence = "A" * 50
    record.flag = flag
    record.reference_id = 0
    record.reference_start = position
    record.mapping_quality = 60
    record.cigar = ((0, 50),)
    record.query_qualities = pysam.qualitystring_to_array("I" * 50)
    if start is not None:
        record.set_tag("st", start, value_type="Z")
    if duration is not None:
        record.set_tag("du", duration, value_type="f")
    return record


@pytest.fixture
def input_bam(tmp_path: Path) -> Path:
    path = tmp_path / "chunk.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 10000}],
        "RG": [
            {
                "ID": "run",
                "DT": "2026-01-01T00:00:00Z",
                "PL": "ONT",
            }
        ],
    }
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        # Completion times are 3, 11, and 21 seconds from @RG DT.
        bam.write(_record("read-1", 100, "2026-01-01T00:00:01Z", 2.0))
        bam.write(_record("read-2", 200, "2026-01-01T00:00:09Z", 2.0))
        # Dorado/minimap2 commonly omit acquisition tags from supplementary hits.
        bam.write(_record("read-2", 300, None, None, flag=2048))
        bam.write(_record("read-3", 400, "2026-01-01T00:00:19Z", 2.0))
    return path


def _names(path: Path) -> list[str]:
    with pysam.AlignmentFile(path, "rb") as bam:
        return [record.query_name for record in bam.fetch(until_eof=True)]


def test_disjoint_completion_windows(input_bam: Path, tmp_path: Path):
    output = tmp_path / "disjoint"
    windows = split_bams(
        SplitConfig(
            inputs=[str(input_bam)],
            output=output,
            interval_ns=10_000_000_000,
            until_ns=20_000_000_000,
            cumulative=False,
            prefix="case",
        )
    )

    assert [window.path.name for window in windows] == [
        "case_000-010sec.bam",
        "case_010-020sec.bam",
    ]
    assert _names(windows[0].path) == ["read-1"]
    assert _names(windows[1].path) == ["read-2", "read-2"]
    assert Path(str(windows[0].path) + ".bai").exists()
    summary = (output / "timeline_summary.tsv").read_text()
    assert "\t1\t1\t50\n" in summary
    assert "\t1\t2\t50\n" in summary


def test_cumulative_windows(input_bam: Path, tmp_path: Path):
    output = tmp_path / "cumulative"
    windows = split_bams(
        SplitConfig(
            inputs=[str(input_bam)],
            output=output,
            interval_ns=10_000_000_000,
            until_ns=20_000_000_000,
            cumulative=True,
            prefix="case",
            index=False,
        )
    )

    assert [window.path.name for window in windows] == [
        "case_010sec.bam",
        "case_020sec.bam",
    ]
    assert _names(windows[0].path) == ["read-1"]
    assert _names(windows[1].path) == ["read-1", "read-2", "read-2"]


def test_existing_outputs_are_protected(input_bam: Path, tmp_path: Path):
    config = SplitConfig(
        inputs=[str(input_bam)],
        output=tmp_path / "protected",
        interval_ns=10_000_000_000,
        until_ns=10_000_000_000,
        cumulative=False,
        index=False,
    )
    split_bams(config)
    with pytest.raises(NanoTimeError, match="output already exists"):
        split_bams(config)

    forced = SplitConfig(**{**config.__dict__, "force": True})
    split_bams(forced)


def test_start_mode_does_not_require_duration(tmp_path: Path):
    path = tmp_path / "no-duration.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 1000}],
        "RG": [{"ID": "run", "DT": "2026-01-01T00:00:00Z"}],
    }
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        bam.write(_record("read", 10, "2026-01-01T00:00:01Z", None))

    windows = split_bams(
        SplitConfig(
            inputs=[str(path)],
            output=tmp_path / "out",
            interval_ns=10_000_000_000,
            until_ns=10_000_000_000,
            cumulative=False,
            timestamp_mode="start",
            index=False,
        )
    )
    assert _names(windows[0].path) == ["read"]


def test_end_mode_requires_duration(tmp_path: Path):
    path = tmp_path / "no-duration.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 1000}],
        "RG": [{"ID": "run", "DT": "2026-01-01T00:00:00Z"}],
    }
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        bam.write(_record("read", 10, "2026-01-01T00:00:01Z", None))

    with pytest.raises(NanoTimeError, match="no du tag"):
        split_bams(
            SplitConfig(
                inputs=[str(path)],
                output=tmp_path / "out",
                interval_ns=10_000_000_000,
                until_ns=10_000_000_000,
                cumulative=False,
            )
        )


def test_arbitrary_cumulative_checkpoints(input_bam: Path, tmp_path: Path):
    windows = split_bams(
        SplitConfig(
            inputs=[str(input_bam)],
            output=tmp_path / "checkpoints",
            interval_ns=None,
            until_ns=None,
            checkpoints_ns=(5_000_000_000, 15_000_000_000, 30_000_000_000),
            cumulative=True,
            prefix="case",
            index=False,
        )
    )
    assert [window.path.name for window in windows] == [
        "case_005sec.bam",
        "case_015sec.bam",
        "case_030sec.bam",
    ]
    assert _names(windows[0].path) == ["read-1"]
    assert _names(windows[1].path) == ["read-1", "read-2", "read-2"]
    assert _names(windows[2].path) == ["read-1", "read-2", "read-2", "read-3"]


def test_yield_checkpoints_include_threshold_crossing_read(input_bam: Path, tmp_path: Path):
    output = tmp_path / "yield"
    windows = split_bams(
        SplitConfig(
            inputs=[str(input_bam)],
            output=output,
            interval_ns=None,
            until_ns=None,
            cumulative=True,
            yield_targets=(50, 100, 150),
            prefix="case",
            index=False,
        )
    )
    assert [window.path.name for window in windows] == [
        "case_50b.bam",
        "case_100b.bam",
        "case_150b.bam",
    ]
    assert _names(windows[0].path) == ["read-1"]
    assert _names(windows[1].path) == ["read-1", "read-2", "read-2"]
    assert _names(windows[2].path) == ["read-1", "read-2", "read-2", "read-3"]
    summary = (output / "timeline_summary.tsv").read_text()
    assert "yield_checkpoint_bases\treached_at_seconds" in summary
    assert "\t50\t3\n" in summary


def test_dry_run_writes_nothing(input_bam: Path, tmp_path: Path):
    output = tmp_path / "does-not-exist"
    config = SplitConfig(
        inputs=[str(input_bam)],
        output=output,
        interval_ns=10_000_000_000,
        until_ns=20_000_000_000,
        cumulative=True,
        index=False,
    )
    execution = dry_run_bams(config)
    assert not output.exists()
    assert execution.read_counts == (1, 2)
    assert execution.record_counts == (1, 3)
    assert execution.base_counts == (50, 100)
    assert execution.estimated_output_bytes > 0

    fast_execution = dry_run_bams(SplitConfig(**{**config.__dict__, "fast": True}))
    assert fast_execution.read_counts == execution.read_counts
    assert fast_execution.record_counts == execution.record_counts


def test_inspect_and_timeline_reports(input_bam: Path):
    inspection = inspect_bams([str(input_bam)])
    assert inspection.scan.unique_reads == 3
    assert inspection.scan.total_records == 4
    assert inspection.scan.reads_with_st == 3
    assert inspection.scan.reads_with_du == 3
    assert inspection.scan.supplementary_without_st == 1
    assert inspection.yields[0] == (10, 3, 150)

    timeline = timeline_bams([str(input_bam)], 10_000_000_000)
    assert [(row.new_reads, row.cumulative_reads, row.new_bases) for row in timeline.rows] == [
        (1, 1, 50),
        (1, 2, 50),
        (1, 3, 50),
    ]


def test_manifest_describes_disjoint_composition(input_bam: Path, tmp_path: Path):
    output = tmp_path / "manifest"
    split_bams(
        SplitConfig(
            inputs=[str(input_bam)],
            output=output,
            interval_ns=10_000_000_000,
            until_ns=20_000_000_000,
            cumulative=False,
            manifest=True,
            index=False,
        )
    )
    manifest = json.loads((output / "timeline_manifest.json").read_text())
    assert manifest["mode"] == "disjoint"
    assert manifest["outputs"][1]["combine_through_checkpoint"] == [
        "sample_000-010sec.bam",
        "sample_010-020sec.bam",
    ]


def test_multiple_acquisition_epochs_are_rejected(tmp_path: Path):
    paths = []
    for index, day in enumerate(("01", "02"), start=1):
        path = tmp_path / f"run-{index}.bam"
        header = {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": 1000}],
            "RG": [{"ID": f"run-{index}", "DT": f"2026-01-{day}T00:00:00Z"}],
        }
        with pysam.AlignmentFile(path, "wb", header=header) as bam:
            bam.write(_record(f"read-{index}", index * 10, f"2026-01-{day}T00:00:01Z", 1.0))
        paths.append(path)

    config = SplitConfig(
        inputs=[str(path) for path in paths],
        output=tmp_path / "out",
        interval_ns=10_000_000_000,
        until_ns=10_000_000_000,
        cumulative=False,
        index=False,
    )
    with pytest.raises(NanoTimeError, match="multiple acquisition epochs"):
        split_bams(config)

    allowed = SplitConfig(**{**config.__dict__, "allow_multiple_runs": True})
    execution = dry_run_bams(allowed)
    assert "input contains 2 acquisition epochs" in execution.warnings


def test_suspicious_clock_warnings(tmp_path: Path):
    path = tmp_path / "clock.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 1000}],
        "RG": [{"ID": "run", "DT": "2026-01-01T00:00:00Z"}],
    }
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        bam.write(_record("early", 10, "2025-12-31T23:59:55Z", 1.0))
        bam.write(_record("late", 20, "2026-01-01T06:00:00Z", 1.0))
    report = inspect_bams([str(path)])
    assert any("before clock zero" in warning for warning in report.scan.warnings)
    assert any("hour gap" in warning for warning in report.scan.warnings)
