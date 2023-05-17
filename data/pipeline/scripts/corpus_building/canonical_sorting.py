"""Methods for canonical sorting."""

from functools import cmp_to_key
from typing import Tuple, Set


node_chunk = Tuple[int, int, Set[int]]


def _canonical_order(node_chunk_a: node_chunk, node_chunk_b: node_chunk) -> int:
    """
    Sort items in canonical sorting order.

    Canonical sorting order is defined according to the order
    of Text-Fabric canonical node order. Read more at:
    https://annotation.github.io/text-fabric/tf/core/nodes.html

    This method is used to create a key sorter.

    :param node_chunk_a: the node chunk to get precendence integer for
    :param node_chunk_b: the node chunk to compare with to get precedence for chunk a
    :return: int representing ascending-order precedence of chunk_a versus chunk_b
    """
    na, prec_a, slotsA = node_chunk_a
    nb, prec_b, slotsB = node_chunk_b

    # compare based on node precedence
    if prec_a > prec_b:
        return -1
    elif prec_b > prec_a:
        return 1

    # compare based on slots
    else:
        # slots are equivalent
        if slotsA == slotsB:
            return 0

        # a is subset of b
        aWithoutB = slotsA - slotsB
        if not aWithoutB:
            return 1

        # b is subset of a
        bWithoutA = slotsB - slotsA
        if not bWithoutA:
            return -1

        # compare based on slots
        aMin = min(aWithoutB)
        bMin = min(bWithoutA)
        return -1 if aMin < bMin else 1


canonical_order = cmp_to_key(_canonical_order)
