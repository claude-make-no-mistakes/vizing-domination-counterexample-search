"""Stream graph6 records through an exact, pure Python Vizing factor filter."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from subprocess import PIPE, Popen


GRAPH6_HEADER = b">>graph6<<"


def decode_graph6(record: bytes) -> tuple[int, ...]:
    """Return the closed neighborhood of each vertex as a bit mask."""
    if not record or min(record) < 63 or max(record) > 126:
        raise ValueError(f"expected a graph6 record, got {record!r}")
    data = bytes(value - 63 for value in record)

    order = data[0]
    order_end = 1
    if order == 63:
        order_start, order_end = (2, 8) if data[1:2] == b"\x3f" else (1, 4)
        if len(data) < order_end:
            raise ValueError(f"truncated graph6 order in {record!r}")
        order = 0
        for digit in data[order_start:order_end]:
            order = (order << 6) | digit

    num_pairs = order * (order - 1) // 2
    if len(data) != order_end + (num_pairs + 5) // 6:
        raise ValueError(f"wrong graph6 payload length for order {order}")
    if num_pairs % 6 and data[-1] & ((1 << (6 - num_pairs % 6)) - 1):
        raise ValueError("graph6 record has nonzero padding")

    closed = [1 << vertex for vertex in range(order)]
    pair = 0
    for high in range(1, order):
        for low in range(high):
            if data[order_end + pair // 6] & (32 >> (pair % 6)):
                closed[high] |= 1 << low
                closed[low] |= 1 << high
            pair += 1
    return tuple(closed)


def connected(closed: tuple[int, ...], universe: int) -> bool:
    if not closed:
        return False
    reached = 0
    pending = 1
    while pending:
        vertex_bit = pending & -pending
        reached |= vertex_bit
        pending = (pending | closed[vertex_bit.bit_length() - 1]) & ~reached
    return reached == universe


def has_dominating_set(closed: tuple[int, ...], universe: int, size: int) -> bool:
    for chosen in combinations(closed, size):
        covered = 0
        for neighborhood in chosen:
            covered |= neighborhood
        if covered == universe:
            return True
    return False


def domination_number(closed: tuple[int, ...], universe: int) -> int:
    for size in range(len(closed)):
        if has_dominating_set(closed, universe, size):
            return size
    return len(closed)


def has_packing(closed: tuple[int, ...], size: int) -> bool:
    for chosen in combinations(closed, size):
        covered = 0
        for neighborhood in chosen:
            if covered & neighborhood:
                break
            covered |= neighborhood
        else:
            return True
    return False


def every_edge_addition_reduces(
    closed: tuple[int, ...], universe: int, smaller_size: int
) -> bool:
    for a, b in combinations(range(len(closed)), 2):
        if closed[a] & (1 << b):
            continue
        new_closed = list(closed)
        new_closed[a] |= 1 << b
        new_closed[b] |= 1 << a
        if not has_dominating_set(tuple(new_closed), universe, smaller_size):
            return False
    return True


def identify_pair(closed: tuple[int, ...], a: int, b: int) -> tuple[int, ...]:
    """Merge vertex b into a and leave bit b unused."""
    a_bit = 1 << a
    b_bit = 1 << b
    new_closed = []
    for vertex, neighborhood in enumerate(closed):
        if vertex == b:
            continue
        if vertex == a:
            neighborhood |= closed[b]
        if neighborhood & b_bit:
            neighborhood = (neighborhood & ~b_bit) | a_bit
        new_closed.append(neighborhood)
    return tuple(new_closed)


def every_identification_reduces(
    closed: tuple[int, ...], universe: int, smaller_size: int
) -> bool:
    for a, b in combinations(range(len(closed)), 2):
        new_universe = universe ^ (1 << b)
        if not has_dominating_set(identify_pair(closed, a, b), new_universe, smaller_size):
            return False
    return True


@dataclass
class Counts:
    records: int = 0
    connected: int = 0
    gamma_at_least_4: int = 0
    packing_at_most_gamma_minus_2: int = 0
    every_edge_addition_reduces: int = 0
    every_identification_reduces: int = 0


def filter_graphs(
    records: Iterable[bytes], counts_by_order: dict[int, Counts]
) -> Iterator[bytes]:
    for line in records:
        record = line.strip().removeprefix(GRAPH6_HEADER)
        if not record:
            continue
        closed = decode_graph6(record)
        counts = counts_by_order.setdefault(len(closed), Counts())
        counts.records += 1
        universe = (1 << len(closed)) - 1
        if not connected(closed, universe):
            continue
        counts.connected += 1
        gamma = domination_number(closed, universe)
        if gamma < 4:
            continue
        counts.gamma_at_least_4 += 1
        # A larger packing contains a packing of size gamma - 1.
        if has_packing(closed, gamma - 1):
            continue
        counts.packing_at_most_gamma_minus_2 += 1
        # If a smaller dominating set exists, it extends to one of size gamma - 1.
        if not every_edge_addition_reduces(closed, universe, gamma - 1):
            continue
        counts.every_edge_addition_reduces += 1
        if not every_identification_reduces(closed, universe, gamma - 1):
            continue
        counts.every_identification_reduces += 1
        yield record


def geng_records(geng: Path, min_order: int, max_order: int) -> Iterator[bytes]:
    for order in range(min_order, max_order + 1):
        with Popen([str(geng.resolve()), "-cq", str(order)], stdout=PIPE) as process:
            assert process.stdout is not None
            yield from process.stdout
            if process.wait():
                raise RuntimeError(f"geng failed at order {order} with status {process.returncode}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--input", type=Path, help="Read graph6 from a file. Defaults to stdin.")
    source.add_argument("--geng", type=Path, help="Stream connected graphs from this nauty geng.")
    parser.add_argument("--min-order", type=int, default=1)
    parser.add_argument("--max-order", type=int)
    args = parser.parse_args()

    counts_by_order: dict[int, Counts] = {}
    if args.geng:
        if args.max_order is None or not 1 <= args.min_order <= args.max_order:
            parser.error("--geng requires 1 <= --min-order <= --max-order")
        records = geng_records(args.geng, args.min_order, args.max_order)
    else:
        if args.max_order is not None or args.min_order != 1:
            parser.error("--min-order and --max-order require --geng")
        records = args.input.open("rb") if args.input else sys.stdin.buffer

    try:
        for record in filter_graphs(records, counts_by_order):
            print(record.decode("ascii"))
    finally:
        if records is not sys.stdin.buffer:
            records.close()

    columns = tuple(Counts.__dataclass_fields__)
    print("order", *columns, sep="\t", file=sys.stderr)
    for order, counts in sorted(counts_by_order.items()):
        print(order, *(getattr(counts, column) for column in columns), sep="\t", file=sys.stderr)


if __name__ == "__main__":
    main()
