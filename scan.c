#include "gtools.h"

#include <stdint.h>

static uint64_t eligible_gamma, eligible_packing, eligible_edges, eligible_pairs;

/* The rows of a closed-neighborhood matrix are sets of nauty vertex bits. */
static int covers(const setword *rows, int n, int next, int left, setword missing)
{
    if (!missing) return 1;
    if (!left) return 0;
    for (int i = next; i <= n - left; ++i)
        if (covers(rows, n, i + 1, left - 1, missing & ~rows[i])) return 1;
    return 0;
}

static int covers_three(const setword *rows, int n, setword all)
{
    for (int i = 0; i < n - 2; ++i) {
        setword missing_i = all & ~rows[i];
        for (int j = i + 1; j < n - 1; ++j) {
            setword missing_ij = missing_i & ~rows[j];
            for (int k = j + 1; k < n; ++k)
                if (!(missing_ij & ~rows[k])) return 1;
        }
    }
    return 0;
}

static int packs(const setword *rows, int n, int next, int left, setword used)
{
    if (!left) return 1;
    for (int i = next; i <= n - left; ++i)
        if (!(used & rows[i]) && packs(rows, n, i + 1, left - 1, used | rows[i]))
            return 1;
    return 0;
}

/* Selecting all future vertices adds at most that many to the domination number. */
int vizing_prune(graph *g, int n, int maxn)
{
    int budget = 3 - (maxn - n);
    if (budget < 1 || n == maxn) return 0;
    setword rows[MAXN];
    for (int i = 0; i < n; ++i) rows[i] = g[i] | bit[i];
    return covers(rows, n, 0, budget, ALLMASK(n));
}

void vizing_output(FILE *out, graph *g, int n)
{
    if (n < 4) return;
    setword rows[MAXN], all = ALLMASK(n);
    for (int i = 0; i < n; ++i) rows[i] = g[i] | bit[i];
    if (covers_three(rows, n, all)) return;
    ++eligible_gamma;
    int gamma = 4;
    while (!covers(rows, n, 0, gamma, all)) ++gamma;
    if (packs(rows, n, 0, gamma - 1, 0)) return;
    ++eligible_packing;

    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            if (!(g[i] & bit[j])) {
                rows[i] |= bit[j];
                rows[j] |= bit[i];
                int lowers = covers(rows, n, 0, gamma - 1, all);
                rows[i] ^= bit[j];
                rows[j] ^= bit[i];
                if (!lowers) return;
            }
    ++eligible_edges;

    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j) {
            setword merged[MAXN];
            const setword remaining = all ^ bit[j];
            for (int k = 0; k < n; ++k)
                merged[k] = (rows[k] & remaining) | ((rows[k] & bit[j]) ? bit[i] : 0);
            merged[i] = (rows[i] | rows[j]) & remaining;
            merged[j] = 0;
            if (!covers(merged, n, 0, gamma - 1, remaining)) return;
        }
    ++eligible_pairs;
    writeg6(out, g, 1, n);
}

void vizing_summary(nauty_counter total, double cpu)
{
    (void) cpu;
    fprintf(stderr, "vizing %llu %llu %llu %llu %llu\n",
            (unsigned long long) total,
            (unsigned long long) eligible_gamma,
            (unsigned long long) eligible_packing,
            (unsigned long long) eligible_edges,
            (unsigned long long) eligible_pairs);
}
