"""Duration and timestamp helpers."""

from __future__ import annotations

import math
import re
from datetime import datetime, timezone

_DURATION_RE = re.compile(
    r"^\s*(?P<number>(?:\d+(?:\.\d*)?|\.\d+))\s*(?P<unit>s|m|h|d)?\s*$",
    re.IGNORECASE,
)
_UNIT_SECONDS = {"s": 1, "m": 60, "h": 3600, "d": 86400}
_SIZE_RE = re.compile(
    r"^\s*(?P<number>(?:\d+(?:\.\d*)?|\.\d+))\s*(?P<unit>[kKmMgGtT]?)\s*(?:[bB])?\s*$"
)
_SIZE_MULTIPLIERS = {"": 1, "k": 1_000, "m": 1_000_000, "g": 1_000_000_000, "t": 1_000_000_000_000}


def parse_duration(value: str) -> float:
    """Parse a positive duration such as 30s, 10m, 1.5h, or 1d."""
    match = _DURATION_RE.match(value)
    if not match:
        raise ValueError(
            f"invalid duration {value!r}; use a number followed by s, m, h, or d"
        )
    number = float(match.group("number"))
    unit = (match.group("unit") or "s").lower()
    seconds = number * _UNIT_SECONDS[unit]
    if seconds <= 0:
        raise ValueError("duration must be greater than zero")
    return seconds


def parse_timestamp(value: str) -> datetime:
    """Parse an ISO-8601 timestamp and normalize it to UTC."""
    text = value.strip()
    if text.endswith(("Z", "z")):
        text = text[:-1] + "+00:00"
    try:
        timestamp = datetime.fromisoformat(text)
    except ValueError as exc:
        raise ValueError(f"invalid ISO-8601 timestamp {value!r}") from exc
    if timestamp.tzinfo is None:
        raise ValueError(
            f"timestamp {value!r} has no UTC offset; use Z or an explicit offset"
        )
    return timestamp.astimezone(timezone.utc)


def to_epoch_ns(timestamp: datetime) -> int:
    """Convert a timezone-aware datetime to integer Unix nanoseconds."""
    return int(round(timestamp.timestamp() * 1_000_000_000))


def format_timestamp(epoch_ns: int) -> str:
    """Format integer Unix nanoseconds as an ISO-8601 UTC timestamp."""
    timestamp = datetime.fromtimestamp(epoch_ns / 1_000_000_000, tz=timezone.utc)
    return timestamp.isoformat(timespec="microseconds").replace("+00:00", "Z")


def duration_ns(value: str) -> int:
    """Parse a duration and return integer nanoseconds."""
    nanoseconds = int(round(parse_duration(value) * 1_000_000_000))
    if nanoseconds <= 0:
        raise ValueError("duration is smaller than one nanosecond")
    return nanoseconds


def parse_size(value: str) -> int:
    """Parse a positive decimal base count such as 100M, 1G, or 2.5Gb."""
    match = _SIZE_RE.match(value)
    if not match:
        raise ValueError(
            f"invalid yield {value!r}; use a number optionally followed by K, M, G, or T"
        )
    bases = float(match.group("number")) * _SIZE_MULTIPLIERS[match.group("unit").lower()]
    if not math.isfinite(bases) or bases <= 0:
        raise ValueError("yield must be greater than zero")
    return int(round(bases))


def format_size(value: int, *, compact: bool = False) -> str:
    """Format a base or byte count using decimal SI units."""
    units = ((1_000_000_000_000, "T"), (1_000_000_000, "G"), (1_000_000, "M"), (1_000, "k"))
    for divisor, suffix in units:
        if value >= divisor:
            number = value / divisor
            text = f"{number:.2f}".rstrip("0").rstrip(".")
            return f"{text}{suffix}{'' if compact else 'b'}"
    return f"{value}{'' if compact else 'b'}"


def format_bytes(value: int) -> str:
    """Format a byte count using decimal SI units."""
    units = (
        (1_000_000_000_000, "TB"),
        (1_000_000_000, "GB"),
        (1_000_000, "MB"),
        (1_000, "kB"),
    )
    for divisor, suffix in units:
        if value >= divisor:
            text = f"{value / divisor:.2f}".rstrip("0").rstrip(".")
            return f"{text}{suffix}"
    return f"{value}B"
