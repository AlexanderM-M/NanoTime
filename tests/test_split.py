from pathlib import Path

import pysam
import pytest

from nanotime.split import NanoTimeError, SplitConfig, split_bams


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
