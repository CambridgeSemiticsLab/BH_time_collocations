"""Module to test canonical sorting."""

import pytest

from corpus_building.canonical_sorting import canonical_order


@pytest.mark.parametrize(
    'original, expected',
    [
        (
            [
                (1, 9, {1, 2, 3}),
                (2, 10, {1, 2, 3}),
                (3, 5, {3, 4}),
                (3, 5, {1, 2}),
                (3, 6, {1, 2, 3}),
                (3, 6, {1, 2})
            ],
            [
                (2, 10, {1, 2, 3}),
                (1, 9, {1, 2, 3}),
                (3, 6, {1, 2, 3}),
                (3, 6, {1, 2}),
                (3, 5, {1, 2}),
                (3, 5, {3, 4})
            ],
        ),
    ],
)
def test__canonical_order(original, expected):
    actual = sorted(original, key=canonical_order)
    assert actual == expected
