"""Static timeline visualization for NanoTime reports."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

from .split import NanoTimeError, TimelineRow
from .timeutil import format_size


def _axis_scale(maximum: int) -> tuple[float, str]:
    for divisor, label in (
        (1_000_000_000_000, "Tb"),
        (1_000_000_000, "Gb"),
        (1_000_000, "Mb"),
        (1_000, "kb"),
    ):
        if maximum >= divisor:
            return float(divisor), label
    return 1.0, "bases"


def _time_scale(seconds: float) -> tuple[float, str]:
    if seconds <= 180:
        return 1.0, "seconds"
    if seconds <= 3 * 60 * 60:
        return 60.0, "minutes"
    return 3600.0, "hours"


def _duration_text(seconds: float) -> str:
    if seconds < 180:
        return f"{seconds:g} seconds"
    if seconds < 3 * 60 * 60:
        return f"{seconds / 60:g} minutes"
    return f"{seconds / 3600:g} hours"


def _bin_text(seconds: float) -> str:
    if seconds < 180:
        return f"{seconds:g}-second"
    if seconds < 3 * 60 * 60:
        return f"{seconds / 60:g}-minute"
    return f"{seconds / 3600:g}-hour"


def save_timeline_plot(
    rows: Sequence[TimelineRow],
    bin_ns: int,
    output: Path,
    *,
    title: str = "NanoTime sequencing timeline",
    dpi: int = 160,
) -> Path:
    """Render cumulative yield and per-bin throughput to a static image."""
    if not rows:
        raise NanoTimeError("timeline contains no bins to plot")
    if output.suffix.lower() not in {".png", ".svg", ".pdf"}:
        raise NanoTimeError("plot output must end in .png, .svg, or .pdf")

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.ticker import FuncFormatter, MaxNLocator
    except ImportError as exc:
        raise NanoTimeError(
            'plotting requires matplotlib; install it with: pip install "nanotime-ont[plot]"'
        ) from exc

    elapsed_seconds = rows[-1].end_ns / 1_000_000_000
    time_divisor, time_label = _time_scale(elapsed_seconds)
    cumulative_divisor, cumulative_label = _axis_scale(rows[-1].cumulative_bases)
    throughput_divisor, throughput_label = _axis_scale(max(row.new_bases for row in rows))

    ends = [row.end_ns / 1_000_000_000 / time_divisor for row in rows]
    centers = [
        ((row.start_ns + row.end_ns) / 2) / 1_000_000_000 / time_divisor
        for row in rows
    ]
    widths = [
        (row.end_ns - row.start_ns) / 1_000_000_000 / time_divisor * 0.82
        for row in rows
    ]
    cumulative = [row.cumulative_bases / cumulative_divisor for row in rows]
    throughput = [row.new_bases / throughput_divisor for row in rows]

    teal = "#007C83"
    purple = "#71569A"
    grid = "#D8DEE4"
    text = "#25313C"
    background = "#FBFCFE"

    figure, (yield_axis, throughput_axis) = plt.subplots(
        2,
        1,
        figsize=(12, 7.2),
        sharex=True,
        gridspec_kw={"height_ratios": (1.45, 1)},
        constrained_layout=True,
    )
    try:
        figure.patch.set_facecolor(background)
        for axis in (yield_axis, throughput_axis):
            axis.set_facecolor(background)
            axis.grid(axis="y", color=grid, linewidth=0.8, alpha=0.8)
            axis.set_axisbelow(True)
            axis.spines["top"].set_visible(False)
            axis.spines["right"].set_visible(False)
            axis.spines["left"].set_color("#AAB4BE")
            axis.spines["bottom"].set_color("#AAB4BE")
            axis.tick_params(colors=text)
            axis.yaxis.set_major_locator(MaxNLocator(nbins=5, min_n_ticks=3))

        line_x = [0.0, *ends]
        line_y = [0.0, *cumulative]
        yield_axis.plot(line_x, line_y, color=teal, linewidth=2.6)
        yield_axis.fill_between(line_x, line_y, color=teal, alpha=0.13)
        yield_axis.scatter(ends[-1], cumulative[-1], color=teal, s=32, zorder=3)
        yield_axis.set_ylabel(f"Cumulative yield ({cumulative_label})", color=text)
        yield_axis.set_ylim(bottom=0)
        yield_axis.set_title("Yield over time", loc="left", fontsize=12, color=text, pad=8)
        yield_axis.annotate(
            format_size(rows[-1].cumulative_bases),
            xy=(ends[-1], cumulative[-1]),
            xytext=(-8, 10),
            textcoords="offset points",
            ha="right",
            color=teal,
            fontsize=10,
            fontweight="bold",
        )

        throughput_axis.bar(
            centers,
            throughput,
            width=widths,
            color=purple,
            edgecolor="none",
            alpha=0.9,
        )
        throughput_axis.set_ylabel(f"New bases ({throughput_label})", color=text)
        throughput_axis.set_xlabel(f"Elapsed time ({time_label})", color=text)
        throughput_axis.set_ylim(bottom=0)
        bin_seconds = bin_ns / 1_000_000_000
        throughput_axis.set_title(
            f"Throughput per {_bin_text(bin_seconds)} bin",
            loc="left",
            fontsize=12,
            color=text,
            pad=8,
        )
        throughput_axis.set_xlim(0, max(ends))
        throughput_axis.xaxis.set_major_locator(MaxNLocator(nbins=9, min_n_ticks=4))
        throughput_axis.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))

        total_reads = rows[-1].cumulative_reads
        summary = (
            f"{format_size(rows[-1].cumulative_bases)} total yield  •  "
            f"{total_reads:,} reads  •  {_duration_text(elapsed_seconds)}"
        )
        figure.suptitle(title, fontsize=18, fontweight="bold", color=text, y=1.035)
        figure.text(0.5, 0.98, summary, ha="center", va="top", fontsize=10.5, color="#5B6772")

        output.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(
            output,
            dpi=dpi,
            bbox_inches="tight",
            facecolor=figure.get_facecolor(),
            metadata={"Title": title, "Creator": "NanoTime"},
        )
    finally:
        plt.close(figure)
    return output
