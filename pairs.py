"""Check fixed graph pairs for a counterexample to Vizing's domination conjecture."""

from __future__ import annotations

import argparse
import json
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import asdict, dataclass, field
from functools import lru_cache
from itertools import combinations
from pathlib import Path
from time import perf_counter
from typing import TYPE_CHECKING, Literal

import networkx as nx

if TYPE_CHECKING:
    from pysat.formula import CNF


Vertex = tuple[int, int]


def _domination_number(neighbors: tuple[frozenset[int], ...]) -> int:
    order = len(neighbors)
    if order == 0:
        return 0

    masks = tuple(
        (1 << u) | sum(1 << v for v in adjacent)
        for u, adjacent in enumerate(neighbors)
    )
    all_vertices = (1 << order) - 1
    largest_neighborhood = max(mask.bit_count() for mask in masks)
    lower_bound = (order + largest_neighborhood - 1) // largest_neighborhood

    for size in range(lower_bound, order + 1):
        for selected in combinations(masks, size):
            covered = 0
            for mask in selected:
                covered |= mask
            if covered == all_vertices:
                return size
    raise AssertionError("Selecting all factor vertices must dominate the factor")


@dataclass(frozen=True)
class Factor:
    """A graph6 record with its adjacency and exact domination number."""

    graph6: str
    neighbors: tuple[frozenset[int], ...] = field(init=False, repr=False)
    gamma: int = field(init=False)

    def __post_init__(self) -> None:
        try:
            graph = nx.from_graph6_bytes(self.graph6.encode("ascii"))
        except (IndexError, ValueError, nx.NetworkXError) as error:
            raise ValueError(f"Invalid graph6 record {self.graph6!r} ({error})") from error
        neighbors = tuple(frozenset(graph[u]) for u in range(len(graph)))
        object.__setattr__(self, "neighbors", neighbors)
        object.__setattr__(self, "gamma", _domination_number(neighbors))


FactorInput = Factor | str | bytes


@lru_cache(maxsize=256)
def _prepare_record(graph6: str) -> Factor:
    return Factor(graph6)


def prepare_factor(graph6: FactorInput) -> Factor:
    """Parse a graph6 record and compute its domination number by enumeration."""
    if isinstance(graph6, Factor):
        return graph6
    if isinstance(graph6, bytes):
        graph6 = graph6.decode("ascii")
    return _prepare_record(graph6.strip().removeprefix(">>graph6<<"))


def dominates_product(
    left: Factor, right: Factor, vertices: Iterable[Sequence[int]]
) -> bool:
    """Check distinct, valid product vertices and domination by direct adjacency."""
    left_order, right_order = len(left.neighbors), len(right.neighbors)
    selected: set[Vertex] = set()
    for vertex in vertices:
        try:
            u, v = vertex
        except (TypeError, ValueError):
            return False
        if type(u) is not int or type(v) is not int:
            return False
        if not (0 <= u < left_order and 0 <= v < right_order) or (u, v) in selected:
            return False
        selected.add((u, v))

    return all(
        any(
            (u == a and (v == b or b in right.neighbors[v]))
            or (v == b and a in left.neighbors[u])
            for a, b in selected
        )
        for u in range(left_order)
        for v in range(right_order)
    )


def verify_counterexample(
    left: FactorInput, right: FactorInput, vertices: Iterable[Sequence[int]]
) -> bool:
    """Verify a dominating set strictly below the product of the factor optima."""
    left_factor, right_factor = prepare_factor(left), prepare_factor(right)
    candidate = tuple(vertices)
    return (
        len(candidate) < left_factor.gamma * right_factor.gamma
        and dominates_product(left_factor, right_factor, candidate)
    )


def product_domination_cnf(left: Factor, right: Factor, bound: int) -> CNF:
    """Encode domination with at most bound vertices in the Cartesian product."""
    from pysat.card import CardEnc, EncType
    from pysat.formula import CNF

    right_order = len(right.neighbors)
    product_order = len(left.neighbors) * right_order
    clauses = CNF()
    if bound < 0:
        clauses.append([])
        return clauses

    for u, left_neighbors in enumerate(left.neighbors):
        for v, right_neighbors in enumerate(right.neighbors):
            vertex = u * right_order + v + 1
            clause = [vertex]
            clause.extend(a * right_order + v + 1 for a in sorted(left_neighbors))
            clause.extend(u * right_order + b + 1 for b in sorted(right_neighbors))
            clauses.append(clause)

    clauses.extend(CardEnc.atmost(
        lits=range(1, product_order + 1), bound=bound, encoding=EncType.seqcounter
    ).clauses)
    return clauses


@dataclass(frozen=True)
class PairResult:
    left_graph6: str
    right_graph6: str
    left_gamma: int
    right_gamma: int
    bound: int
    status: Literal["counterexample", "unsat", "unknown"]
    witness: tuple[Vertex, ...] | None
    solver: str
    conflict_budget: int | None
    conflicts: int | None
    solve_seconds: float
    proof_lines: int | None


def check_pair(
    left: FactorInput,
    right: FactorInput,
    *,
    solver: str | None = None,
    conflict_budget: int | None = None,
    cnf_path: str | Path | None = None,
    proof_path: str | Path | None = None,
) -> PairResult:
    """Check one pair. UNSAT is a solver verdict and proof output needs replay.

    The default solver is Glucose4 for complete calls and CaDiCaL 1.9.5 for
    limited calls. The conflict budget depends on the backend's implementation.
    """
    from pysat.solvers import Solver

    if conflict_budget is not None and conflict_budget <= 0:
        raise ValueError("The SAT conflict budget must be positive")
    if solver is None:
        solver = "g4" if conflict_budget is None else "cadical195"

    left_factor, right_factor = prepare_factor(left), prepare_factor(right)
    right_order = len(right_factor.neighbors)
    product_order = len(left_factor.neighbors) * right_order
    bound = left_factor.gamma * right_factor.gamma - 1
    formula = product_domination_cnf(left_factor, right_factor, bound)
    if cnf_path is not None:
        formula.to_file(str(cnf_path))

    witness = None
    proof_lines = None
    status: Literal["counterexample", "unsat", "unknown"]
    with Solver(name=solver, bootstrap_with=formula.clauses, with_proof=proof_path is not None) as sat:
        started = perf_counter()
        if conflict_budget is None:
            solved = sat.solve()
        else:
            sat.conf_budget(conflict_budget)
            solved = sat.solve_limited()
        solve_seconds = perf_counter() - started
        stats = sat.accum_stats()
        # Drain the proof stream before closing the backend, including after an
        # interrupted solve. PySAT's Windows CaDiCaL builds can crash otherwise.
        proof = sat.get_proof() if proof_path is not None else None
        if solved is True:
            witness = tuple(sorted(
                divmod(lit - 1, right_order) for lit in sat.get_model()
                if 1 <= lit <= product_order
            ))
            if not verify_counterexample(left_factor, right_factor, witness):
                raise RuntimeError("The SAT model failed direct counterexample verification")
            status = "counterexample"
        elif solved is False:
            status = "unsat"
            if proof_path is not None:
                if proof is None:
                    raise RuntimeError(f"Solver {solver!r} did not return an UNSAT proof")
                Path(proof_path).write_text("".join(f"{line}\n" for line in proof), encoding="ascii")
                proof_lines = len(proof)
        else:
            status = "unknown"

    return PairResult(
        left_graph6=left_factor.graph6,
        right_graph6=right_factor.graph6,
        left_gamma=left_factor.gamma,
        right_gamma=right_factor.gamma,
        bound=bound,
        status=status,
        witness=witness,
        solver=solver,
        conflict_budget=conflict_budget,
        conflicts=None if stats is None else stats.get("conflicts"),
        solve_seconds=solve_seconds,
        proof_lines=proof_lines,
    )


def check_pairs(
    pairs: Iterable[tuple[FactorInput, FactorInput]],
    *,
    solver: str | None = None,
    conflict_budget: int | None = None,
) -> Iterator[PairResult]:
    """Yield results in input order, reusing recently prepared factors."""
    for left, right in pairs:
        yield check_pair(left, right, solver=solver, conflict_budget=conflict_budget)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("left", help="First graph6 record")
    parser.add_argument("right", help="Second graph6 record")
    parser.add_argument("--solver", help="PySAT solver name (default g4 or cadical195 with a budget)")
    parser.add_argument("--conflicts", type=int, help="Request a conflict budget from the SAT backend")
    parser.add_argument("--cnf", type=Path, help="Write the DIMACS formula to this file")
    parser.add_argument("--proof", type=Path, help="Write an UNSAT proof log to this file")
    parser.add_argument("--verify", type=Path, help="Verify the product vertices in a JSON witness file")
    args = parser.parse_args()

    try:
        if args.verify is not None:
            data = json.loads(args.verify.read_text(encoding="utf-8"))
            witness = data.get("witness") if isinstance(data, dict) else data
            valid = isinstance(witness, list) and verify_counterexample(args.left, args.right, witness)
            print(json.dumps({"valid_counterexample": valid}))
            return 0 if valid else 1

        result = check_pair(
            args.left, args.right, solver=args.solver, conflict_budget=args.conflicts,
            cnf_path=args.cnf, proof_path=args.proof,
        )
    except (ValueError, OSError) as error:
        parser.error(str(error))
    print(json.dumps(asdict(result)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
