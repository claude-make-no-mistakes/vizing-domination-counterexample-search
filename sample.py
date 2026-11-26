"""Generate domination-critical candidates from randomly saturated cycles."""

from argparse import ArgumentParser
from itertools import combinations, combinations_with_replacement
from pathlib import Path
from random import Random

import networkx as nx

from independent import decode_graph6, has_dominating_set, has_packing


def neighborhoods(graph):
    record = nx.to_graph6_bytes(graph, header=False).strip()
    return record, decode_graph6(record)


def critical_candidate(order, gamma, random):
    # Adding one edge can reduce domination number by at most one.
    graph = nx.cycle_graph(order)
    while True:
        edges = list(nx.non_edges(graph))
        random.shuffle(edges)
        for a, b in edges:
            graph.add_edge(a, b)
            _, rows = neighborhoods(graph)
            if has_dominating_set(rows, (1 << len(rows)) - 1, gamma - 1):
                graph.remove_edge(a, b)

        pairs = list(combinations(graph, 2))
        random.shuffle(pairs)
        for a, b in pairs:
            smaller = nx.contracted_nodes(graph, a, b, self_loops=False)
            _, rows = neighborhoods(smaller)
            if not has_dominating_set(rows, (1 << len(rows)) - 1, gamma - 1):
                graph = nx.convert_node_labels_to_integers(smaller)
                break
        else:
            record, rows = neighborhoods(graph)
            if has_packing(rows, gamma - 1):
                return None
            return record


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=20261126)
    parser.add_argument("--gamma", type=int, default=5)
    parser.add_argument("--start", type=int, default=13)
    parser.add_argument("--stop", type=int, default=16)
    parser.add_argument("--trials", type=int, default=8)
    parser.add_argument("--check", type=Path, help="check all pairs from a saved sample")
    parser.add_argument("--include", action="append", default=[], help="add a factor to --check")
    parser.add_argument("--conflicts", type=int, default=50_000)
    args = parser.parse_args()
    if args.check is not None:
        from dataclasses import asdict
        import json

        from pairs import check_pairs

        records = args.check.read_text(encoding="ascii").splitlines() + args.include
        for result in check_pairs(combinations_with_replacement(records, 2),
                                  conflict_budget=args.conflicts):
            print(json.dumps(asdict(result)), flush=True)
        return
    if not (4 <= args.gamma <= (args.start + 2) // 3 <= (args.stop + 2) // 3
            and args.start <= args.stop <= 62 and args.trials >= 1):
        parser.error("require 4 <= gamma <= ceil(start/3), start <= stop <= 62, and trials >= 1")
    random = Random(args.seed)
    records = set()
    for order in range(args.start, args.stop + 1):
        for _ in range(args.trials):
            result = critical_candidate(order, args.gamma, random)
            if result is not None:
                records.add(result)
    for record in sorted(records):
        print(record.decode("ascii"))


if __name__ == "__main__":
    main()
