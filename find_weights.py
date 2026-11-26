"""Find rational weighted-block certificates for small graph6 factors."""

from argparse import ArgumentParser
from fractions import Fraction
import json
from math import lcm
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

from independent import domination_number
from reception import closed_neighborhoods, find_reception, subsets, verify_weights


def find_weights(rows, gamma):
    blocks = np.arange(1, 1 << len(rows))
    choices = [
        (selected, covered, size)
        for size in range(1, gamma)
        for selected, covered in subsets(rows, size)
    ]
    active = [((blocks & selected) != 0) | ((blocks & covered) == blocks)
              for selected, covered, _ in choices]
    active.append(np.ones(len(blocks), dtype=bool))
    result = linprog(
        -np.ones(len(blocks)), A_ub=active,
        b_ub=[size for _, _, size in choices] + [gamma],
        method="highs",
    )
    if not result.success:
        raise RuntimeError(result.message)

    # The linear program only proposes weights. Reject any inexact candidate.
    approximations = [(int(block), Fraction(float(value)).limit_denominator(10_000))
                      for block, value in zip(blocks, result.x) if value > 1e-8]
    denominator = lcm(*(weight.denominator for _, weight in approximations))
    numerators = [(block, int(weight * denominator)) for block, weight in approximations
                  if weight]
    verify_weights(rows, denominator, numerators)
    return denominator, numerators


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("files", type=Path, nargs="+")
    parser.add_argument("--try-fair-first", action="store_true")
    args = parser.parse_args()
    certificates = {}
    for source in args.files:
        for record in source.read_text(encoding="ascii").splitlines():
            rows = closed_neighborhoods(record)
            gamma = domination_number(tuple(rows), (1 << len(rows)) - 1)
            reception = find_reception(rows, gamma) if args.try_fair_first else None
            if reception is not None:
                denominator, weights = 1, [(block, 1) for block in reception]
            else:
                denominator, weights = find_weights(rows, gamma)
            certificates[record] = {"gamma": gamma, "denominator": denominator,
                                    "weights": weights}
    print(json.dumps(certificates, separators=(",", ":"), ensure_ascii=True))


if __name__ == "__main__":
    main()
