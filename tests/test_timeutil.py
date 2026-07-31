from datetime import timezone

import pytest

from nanotime.timeutil import duration_ns, parse_duration, parse_timestamp


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
