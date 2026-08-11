from datetime import timezone

import pytest

from nanotime.timeutil import (
    duration_ns,
    format_bytes,
    format_size,
    parse_duration,
    parse_size,
    parse_timestamp,
)


@pytest.mark.parametrize(
    ("text", "seconds"),
    [
        ("30", 30),
        ("30s", 30),
        ("10m", 600),
        ("1.5h", 5400),
        ("1d", 86400),
    ],
)
def test_parse_duration(text, seconds):
    assert parse_duration(text) == seconds


@pytest.mark.parametrize("text", ["", "ten minutes", "10x", "0m"])
def test_reject_invalid_duration(text):
    with pytest.raises(ValueError):
        parse_duration(text)


def test_parse_timestamp_normalizes_to_utc():
    timestamp = parse_timestamp("2026-07-29T13:43:45.335073+02:00")
    assert timestamp.tzinfo == timezone.utc
    assert timestamp.isoformat() == "2026-07-29T11:43:45.335073+00:00"


def test_duration_ns():
    assert duration_ns("0.5s") == 500_000_000


@pytest.mark.parametrize(
    ("text", "bases"),
    [("100M", 100_000_000), ("1Gb", 1_000_000_000), ("2.5k", 2_500)],
)
def test_parse_size(text, bases):
    assert parse_size(text) == bases


def test_format_size():
    assert format_size(1_500_000_000) == "1.5Gb"
    assert format_size(100_000_000, compact=True) == "100M"
    assert format_bytes(1_500_000_000) == "1.5GB"
