"""Replay the factor ledger and check the coverage of the recorded pair proofs."""

from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction
from itertools import combinations_with_replacement
import json
from pathlib import Path
from subprocess import run

from independent import decode_graph6, domination_number, filter_graphs
from reception import verify_weights


ROOT = Path(__file__).resolve().parent


def load_json(name):
    return json.loads((ROOT / name).read_text(encoding="ascii"))


def load_jsonl(name):
    return [json.loads(line) for line in (ROOT / name).read_text(encoding="ascii").splitlines()]


def factors_for(order):
    path = ROOT / f"candidates-{order}.g6"
    return path.read_bytes().splitlines() if path.is_file() else []


def replay_scan(executable, mode, table):
    for key, expected in table[mode].items():
        order = int(key)
        cmd = [str(executable), "-cq"]
        if mode == "prefix_degree":
            cmd.append(f"-D{order - 5}")
        result = run(cmd + [str(order)], capture_output=True, check=True)
        summary = next(line for line in result.stderr.splitlines()
                       if line.startswith(b"vizing "))
        actual = list(map(int, summary.split()[1:]))
        if actual != expected or sorted(result.stdout.splitlines()) != sorted(factors_for(order)):
            raise RuntimeError(f"scanner mismatch at order {order}")
        print(f"scanner confirmed order {order}", flush=True)


def verify_ledger(proof_log):
    table = load_json("scan-counts.json")
    factors = []
    for order, counts in table["unpruned"].items():
        retained = factors_for(int(order))
        if len(set(retained)) != counts[-1] or list(filter_graphs(retained, {})) != retained:
            raise RuntimeError(f"critical-factor replay failed at order {order}")
        factors.extend(record.decode("ascii") for record in retained)

    certificates = load_json("weights-lp.json")
    if set(certificates) != set(factors):
        raise RuntimeError("weight records do not match the critical-factor ledger")
    bounds = {}
    for record, certificate in certificates.items():
        rows = decode_graph6(record.encode("ascii"))
        gamma = domination_number(rows, (1 << len(rows)) - 1)
        if gamma != certificate["gamma"]:
            raise RuntimeError(f"wrong domination number for {record}")
        bounds[record] = (gamma, Fraction(*verify_weights(
            rows, certificate["denominator"], certificate["weights"])))
    residual = []
    for a, b in combinations_with_replacement(certificates, 2):
        k, left_weight = bounds[a]
        l, right_weight = bounds[b]
        if left_weight * l <= k * l - 1 and right_weight * k <= k * l - 1:
            residual.append((a, b))

    proof_rows = [json.loads(line) for line in proof_log.read_text(encoding="ascii").splitlines()]
    if len(proof_rows) != len(residual):
        raise RuntimeError("wrong number of residual proof records")
    for index, (pair, row) in enumerate(zip(residual, proof_rows)):
        if ((row["left"], row["right"]) != pair or row["index"] != index
                or row["status"] != "verified" or row["proof_lines"] < 1):
            raise RuntimeError(f"incorrect proof record for residual pair {index}")
    print(f"{len(bounds)} exact factor certificates and {len(residual)} residual proof records")
    for order in (10, 11):
        smallest = min(weight / gamma for record, (gamma, weight) in bounds.items()
                       if len(decode_graph6(record.encode("ascii"))) == order)
        print(f"order {order} minimum certified coefficient {smallest}")

    sample = load_jsonl("sample-pairs.jsonl")
    sample_factors = (ROOT / "sample-gamma5.g6").read_text(encoding="ascii").splitlines()
    exceptions = [g for g in bounds if len(decode_graph6(g.encode("ascii"))) == 10
                  and bounds[g][1] < 4]
    expected = list(combinations_with_replacement(sample_factors + exceptions, 2))
    if [(row["left_graph6"], row["right_graph6"]) for row in sample] != expected:
        raise RuntimeError("sample ledger does not cover the specified pair list")
    statuses = Counter(row["status"] for row in sample)
    print(f"{len(sample)} sampled pairs, {dict(statuses)}")
    return table


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--proof-log", type=Path, default=ROOT / "residual-pairs.jsonl")
    parser.add_argument("--scan", type=Path, help="compiled geng with the scanner hooks")
    parser.add_argument("--mode", choices=["unpruned", "prefix", "prefix_degree"],
                        default="unpruned")
    args = parser.parse_args()
    table = verify_ledger(args.proof_log)
    if args.scan is not None:
        replay_scan(args.scan.resolve(strict=True), args.mode, table)


if __name__ == "__main__":
    main()
