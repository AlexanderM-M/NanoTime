"""Regenerate the README timeline preview after installing NanoTime[plot]."""

from __future__ import annotations

import math
from pathlib import Path

from nanotime.plot import save_timeline_plot
from nanotime.split import MINUTE_NS, TimelineRow


def example_rows() -> list[TimelineRow]:
    rows: list[TimelineRow] = []
    cumulative_reads = 0
    cumulative_bases = 0
    bin_ns = 5 * MINUTE_NS
    for index in range(24):
        trend = 205 + 5.5 * index
        variation = 38 * math.sin(index * 0.72) + 19 * math.cos(index * 0.31)
        pause = 0.32 if index in {13, 14} else 1.0
        new_bases = int(max(35, (trend + variation) * pause) * 1_000_000)
        new_reads = round(new_bases / 11_800)
        cumulative_bases += new_bases
        cumulative_reads += new_reads
        rows.append(
            TimelineRow(
                index=index,
                start_ns=index * bin_ns,
                end_ns=(index + 1) * bin_ns,
                new_reads=new_reads,
                cumulative_reads=cumulative_reads,
                new_bases=new_bases,
                cumulative_bases=cumulative_bases,
            )
        )
    return rows


if __name__ == "__main__":
    destination = Path(__file__).parent / "images" / "timeline-example.png"
    save_timeline_plot(example_rows(), 5 * MINUTE_NS, destination, dpi=170)
    print(destination)
