from pathlib import Path

import pysam

from nanotime.cli import main


def _bam(path: Path) -> Path:
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 1000}],
        "RG": [{"ID": "run", "DT": "2026-01-01T00:00:00Z"}],
    }
    record = pysam.AlignedSegment()
    record.query_name = "read"
    record.query_sequence = "A" * 100
    record.flag = 0
    record.reference_id = 0
    record.reference_start = 10
    record.mapping_quality = 60
    record.cigar = ((0, 100),)
    record.query_qualities = pysam.qualitystring_to_array("I" * 100)
    record.set_tag("st", "2026-01-01T00:00:01Z", value_type="Z")
    record.set_tag("du", 1.0, value_type="f")
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        bam.write(record)
    return path


def test_cli_inspect_and_timeline(tmp_path: Path, capsys):
    bam = _bam(tmp_path / "input.bam")
    assert main(["inspect", str(bam), "--no-progress"]) == 0
    inspected = capsys.readouterr().out
    assert "NanoTime timeline inspection" in inspected
    assert "Unique reads            1" in inspected
    assert "Suggested command:" in inspected

    timeline = tmp_path / "timeline.tsv"
    assert main(
        ["timeline", str(bam), "--bin", "1m", "--output", str(timeline), "--no-progress"]
    ) == 0
    assert timeline.read_text().splitlines() == [
        "bin_start_seconds\tbin_end_seconds\tnew_reads\tcumulative_reads\tnew_bases\tcumulative_bases",
        "0\t60\t1\t1\t100\t100",
    ]

    image = tmp_path / "timeline.png"
    assert main(
        ["plot", str(bam), "--bin", "1m", "--output", str(image), "--no-progress"]
    ) == 0
    assert image.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    assert image.stat().st_size > 10_000


def test_cli_plot_rejects_unknown_image_format(tmp_path: Path, capsys):
    bam = _bam(tmp_path / "input.bam")
    assert main(
        ["plot", str(bam), "--bin", "1m", "--output", str(tmp_path / "plot.jpg"), "--no-progress"]
    ) == 2
    assert "must end in .png, .svg, or .pdf" in capsys.readouterr().err


def test_cli_dry_run_creates_no_output(tmp_path: Path, capsys):
    bam = _bam(tmp_path / "input.bam")
    output = tmp_path / "not-created"
    assert main(
        [
            "split",
            str(bam),
            "--yield",
            "100",
            "--output",
            str(output),
            "--dry-run",
            "--no-progress",
        ]
    ) == 0
    captured = capsys.readouterr()
    assert "sample_100b.bam" in captured.out
    assert "checkpoint\treached_at" in captured.out
    assert "Estimated total BAM size" in captured.out
    assert not output.exists()
