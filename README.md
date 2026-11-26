# Vizing domination counterexample search

Vizing's conjecture asks whether the domination number of the Cartesian
product of any finite simple graphs is at least the product of their domination
numbers. This search found no counterexample.

- The inequality holds if both factors have at most 11 vertices.
- If one factor has at most 9 vertices, the inequality holds for an arbitrary
  partner. The exact inequality also holds if one factor has at most 10
  vertices and the other has domination number at most 4.
- For a factor of order at most 11 and an arbitrary partner, the product
  domination number is at least $61/66$ of the conjectured product.
  The factor is $19/20$ for order at most 10.

[Mathematical note](note.pdf) · [LaTeX source](note.tex)

The connected order-11 scan visited 1,006,700,565 unlabeled graphs. Published
reductions leave 63 candidates through order 11, each with domination number
4 and packing number 2. Rational weighted certificates rule out 1,971 of
their 2,016 unordered pairs. The remaining 45 pairs have DRAT refutations
checked with DRAT-trim. Independent Python enumeration agrees through order 10
and rechecks all 58 retained order-11 factors.

## Replay

Use Python 3.10 or later and the dependencies in [requirements.txt](requirements.txt)
for the SAT tools. The weight checker and independent graph6 filter use the
standard library. The optional weight search uses NumPy and SciPy. A replay
checks the recorded rational weights directly, so it can use the supplied
weights without a numerical optimizer.

Run the lightweight ledger audit and check for agreement with the retained
factors, rational inequalities and proof-log coverage.

~~~sh
python audit.py
~~~

Build the scanner with the source distribution of [nauty 2.9.3](https://pallini.di.uniroma1.it/).
Set NAUTY to the unpacked directory after running nauty's configuration step.
The following command uses a POSIX shell with GCC. Use a C compiler with
32-bit setwords or a GCC-compatible toolchain on Windows and give the
binary an .exe suffix there.

~~~sh
cc -O3 -DMAXN=16 -DWORDSIZE=32 -DOUTPROC=vizing_output -DSUMMARY=vizing_summary \
  -I "$NAUTY" "$NAUTY/geng.c" "$NAUTY/nauty.c" "$NAUTY/nautil.c" \
  "$NAUTY/naugraph.c" "$NAUTY/schreier.c" "$NAUTY/naurng.c" "$NAUTY/gtools.c" \
  scan.c -o scan
python audit.py --scan ./scan
~~~

The five counters in [scan-counts.json](scan-counts.json) are generated graphs,
domination number at least 4, packing number at most domination number minus
2, all missing-edge additions lowering domination, and all unordered
identifications lowering domination. Add `-DPREPRUNE=vizing_prune` to the
build for the faster variant. Replay that binary with
`python audit.py --scan ./scan --mode prefix_degree`. The independent full
order-10 run uses an ordinary geng executable without the hooks.

~~~sh
python independent.py --geng /path/to/geng --max-order 10
python residual.py weights-lp.json --drat-checker /path/to/drat-trim > replay.jsonl
~~~

The Boolean formulas and proofs are generated during replay. The supplied
[pair ledger](residual-pairs.jsonl) records the checker verdict and proof size
for each pair. Replay uses disk space for one proof at a time and streams the
text trace to [DRAT-trim](https://github.com/marijnheule/drat-trim).
Without `--drat-checker` the command reports the SAT solver's verdict.
Use `python audit.py --proof-log replay.jsonl` to check that a regenerated
ledger covers every residual pair. The optional `--canonical` mode adds
the full-projection condition and lexicographic constraints for verified
product automorphisms. Its ledger rows identify that encoding.

## Larger-factor search

The seeded edge-addition and identification search found 22 factors of orders
13 through 16 with domination number 5 and packing number 3. The saved
[sample](sample-gamma5.g6) was tested together with the order-10 graph with
weighted coefficient $19/5$. Of the 276 bounded checks, five returned
UNSAT and 271 were inconclusive at 50,000 conflicts. The sample verdicts
were recorded separately from the 45 proof-checked claims.

~~~sh
python sample.py > regenerated-sample.g6
exception=$(head -n 1 candidates-10.g6)
python sample.py --check sample-gamma5.g6 --include "$exception"
~~~

`pairs.py` also takes two graph6 records and optionally exports a DIMACS
formula and an UNSAT proof. Its witness verifier recomputes both domination
numbers and checks product domination directly.

## Prior art

The September 2026 revision by Aliabadi and Krop proves a universal $2/3$
bound. Indexed web and direct X searches on the ordinary domination conjecture
located statements and partial results. The X audit cannot cover unindexed
posts or private work. [September preprint](https://arxiv.org/abs/2607.01109v3)

The indexed searches also located computational sum-of-squares certificates
for particular graph classes. I found no prior complete verification of the
order-11 pair range in the sources inspected. [Optimization-based proofs](https://arxiv.org/abs/2003.04021) · [Gröbner-basis proofs](https://arxiv.org/abs/2112.04007)

An August 2026 Zenodo deposit claims a Lean resolution. Its final theorem
assumes the desired product lower bound. Its asserted fractional product
identity also fails for $K_2\square K_2=C_4$, whose fractional domination
number is $4/3$ while both factors have fractional domination number 1.
[Deposit and Lean source](https://zenodo.org/records/21854681)
