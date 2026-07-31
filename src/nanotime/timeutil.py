"""Duration and timestamp helpers."""

from __future__ import annotations

import re
from datetime import datetime, timezone

_DURATION_RE = re.compile(
    r"^\s*(?P<number>(?:\d+(?:\.\d*)?|\.\d+))\s*(?P<unit>s|m|h|d)?\s*$",
    re.IGNORECASE,
)
_UNIT_SECONDS = {"s": 1, "m": 60, "h": 3600, "d": 86400}


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
    return int(round(parse_duration(value) * 1_000_000_000))
