"""Generate and verify SAT certificates for pairs left by the weighted bound."""

from argparse import ArgumentParser
from itertools import combinations, combinations_with_replacement
import json
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory

import networkx as nx
from pysat.card import CardEnc, EncType
from pysat.solvers import Solver

from pairs import prepare_factor, product_domination_cnf, verify_counterexample
from reception import closed_neighborhoods, subsets, verify_weights


def remaining_pairs(certificates):
    cases = []
    for graph, certificate in certificates.items():
        rows = closed_neighborhoods(graph)
        weight, denominator = verify_weights(
            rows, certificate["denominator"], certificate["weights"]
        )
        gamma = prepare_factor(graph).gamma
        if gamma != certificate["gamma"]:
            raise ValueError(f"incorrect factor domination number for {graph}")
        cases.append((graph, rows, gamma, weight, denominator, certificate["weights"]))
    for first, second in combinations_with_replacement(cases, 2):
        bound = first[2] * second[2] - 1
        if first[3] * second[2] <= bound * first[4] and second[3] * first[2] <= bound * second[4]:
            yield first, second


def add_fiber_constraints(cnf, factor, partner, variable):
    _, rows, gamma, total, denominator, weights = factor
    _, partner_rows, other_gamma, _, _, _ = partner
    bound = gamma * other_gamma - 1
    slack = denominator * bound - total * other_gamma
    max_size = (total + slack) // denominator
    eligible_blocks = [(block, other_gamma + slack // weight) for block, weight in weights
                       if other_gamma + slack // weight < len(partner_rows)]
    activated = [[] for _ in eligible_blocks]
    for h in range(len(partner_rows)):
        vertices = [variable(g, h) for g in range(len(rows))]
        cnf.extend([-vertex for vertex in choice]
                   for choice in combinations(vertices, max_size + 1))
        flags = [cnf.nv + index + 1 for index in range(len(eligible_blocks))]
        for index, flag in enumerate(flags):
            activated[index].append(flag)
        for size in range(max_size + 1):
            for selected, covered in subsets(rows, size):
                active_blocks = {block for block, _ in weights
                                 if (block & selected) or (block & covered) == block}
                cost = sum(weight for block, weight in weights if block in active_blocks)
                negate_type = [-vertex if (selected >> g) & 1 else vertex
                               for g, vertex in enumerate(vertices)]
                if size * denominator - cost > slack:
                    cnf.append(negate_type)
                else:
                    for flag, (block, _) in zip(flags, eligible_blocks):
                        if block in active_blocks:
                            cnf.append(negate_type + [flag])
        # Reserve flag IDs even if some flags never enter a clause.
        if flags:
            cnf.nv = max(cnf.nv, flags[-1])
    for flags, (_, limit) in zip(activated, eligible_blocks):
        cnf.extend(CardEnc.atmost(lits=flags, bound=limit, top_id=cnf.nv,
                                encoding=EncType.seqcounter).clauses)


def add_lex_leader(cnf, permutation):
    moved = [(old, new) for old, new in enumerate(permutation, start=1) if old != new]
    on_equal_prefix = []
    for index, (left, right) in enumerate(moved):
        cnf.append(on_equal_prefix + [-left, right])
        if index + 1 < len(moved):
            follows = cnf.nv + 1
            cnf.append(on_equal_prefix + [left, right, follows])
            cnf.append(on_equal_prefix + [-left, -right, follows])
            on_equal_prefix = [-follows]


def factor_permutations(factor):
    record, rows, _, _, _, _ = factor
    graph = nx.from_graph6_bytes(record.encode("ascii"))
    for mapping in nx.algorithms.isomorphism.GraphMatcher(graph, graph).isomorphisms_iter():
        permutation = tuple(mapping[vertex] for vertex in range(len(rows)))
        if sorted(permutation) != list(range(len(rows))) or any(
            sum(1 << permutation[j] for j in range(len(rows)) if row & (1 << j))
            != rows[permutation[i]] for i, row in enumerate(rows)
        ):
            raise ValueError("invalid factor automorphism")
        yield permutation


def add_product_symmetry(cnf, first, second):
    height, width = len(first[1]), len(second[1])
    for left in factor_permutations(first):
        for right in factor_permutations(second):
            add_lex_leader(cnf, [left[i] * width + right[j] + 1
                                 for i in range(height) for j in range(width)])
            if first[0] == second[0]:
                add_lex_leader(cnf, [right[j] * width + left[i] + 1
                                     for i in range(height) for j in range(width)])


def add_full_projection(cnf, height, width):
    # An empty row and an empty column would leave their intersection undominated.
    use_rows = cnf.nv + 1
    for i in range(height):
        cnf.append([-use_rows] + [i * width + j + 1 for j in range(width)])
    for j in range(width):
        cnf.append([use_rows] + [i * width + j + 1 for i in range(height)])


def solve_pair(first, second, drat_checker, canonical=False):
    left, right = prepare_factor(first[0]), prepare_factor(second[0])
    width = len(second[1])
    cnf = product_domination_cnf(left, right, left.gamma * right.gamma - 1)
    add_fiber_constraints(cnf, first, second, lambda g, h: g * width + h + 1)
    add_fiber_constraints(cnf, second, first, lambda h, g: g * width + h + 1)
    if canonical:
        add_full_projection(cnf, len(first[1]), width)
        add_product_symmetry(cnf, first, second)
    with Solver(name="cadical195", bootstrap_with=cnf.clauses,
                with_proof=drat_checker is not None) as solver:
        if solver.solve():
            witness = [divmod(lit - 1, width) for lit in solver.get_model()
                       if 1 <= lit <= len(first[1]) * width]
            if not verify_counterexample(left, right, witness):
                raise RuntimeError("invalid counterexample decoded from CNF")
            return {"status": "counterexample", "witness": witness}
        conflicts = solver.accum_stats()["conflicts"]
        if drat_checker is None:
            return {"status": "unsat", "conflicts": conflicts}
        proof = solver.get_proof()
    with TemporaryDirectory(prefix="vizing-drat-") as directory:
        path = Path(directory)
        form = path / "pair.cnf"
        cnf.to_file(str(form))
        checked = run([str(drat_checker), str(form)], input="\n".join(proof) + "\n",
                      capture_output=True, text=True)
        if "s VERIFIED" not in checked.stdout.splitlines():
            raise RuntimeError(checked.stdout + checked.stderr)
        return {"status": "verified", "conflicts": conflicts, "proof_lines": len(proof)}


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("certificates", type=Path)
    parser.add_argument("--drat-checker", type=Path)
    parser.add_argument("--first", type=int, default=0, help="zero-based first pair index")
    parser.add_argument("--stop", type=int, help="exclusive pair index")
    parser.add_argument("--canonical", action="store_true",
                        help="apply product automorphism lex constraints")
    args = parser.parse_args()
    checker = args.drat_checker.resolve(strict=True) if args.drat_checker else None
    cases = remaining_pairs(json.loads(args.certificates.read_text(encoding="ascii")))
    for index, (first, second) in enumerate(cases):
        if index < args.first:
            continue
        if args.stop is not None and index >= args.stop:
            break
        result = solve_pair(first, second, checker, args.canonical)
        if args.canonical:
            result["encoding"] = "canonical"
        print(json.dumps({"index": index, "left": first[0], "right": second[0], **result}),
              flush=True)


if __name__ == "__main__":
    main()
