"""Search for and check fair receptions in small graph6 factors."""

from itertools import combinations


def closed_neighborhoods(record):
    if not record:
        raise ValueError("empty graph6 record")
    n = ord(record[0]) - 63
    if not 1 <= n <= 62 or len(record) != 1 + (n * (n - 1) // 2 + 5) // 6:
        raise ValueError("expected a graph6 record with order 1 through 62")
    payload = [ord(char) - 63 for char in record[1:]]
    if any(not 0 <= digit < 64 for digit in payload):
        raise ValueError("invalid graph6 character")
    rows = [1 << vertex for vertex in range(n)]
    for offset, (right, left) in enumerate(
        (right, left) for right in range(n) for left in range(right)
    ):
        if payload[offset // 6] & (1 << (5 - offset % 6)):
            rows[left] |= 1 << right
            rows[right] |= 1 << left
    return rows


def subsets(rows, size):
    for chosen in combinations(range(len(rows)), size):
        selected = covered = 0
        for vertex in chosen:
            selected |= 1 << vertex
            covered |= rows[vertex]
        yield selected, covered


def verify_reception(rows, blocks):
    full = (1 << len(rows)) - 1
    used = 0
    for block in blocks:
        if block <= 0 or block & (used | ~full):
            return False
        used |= block
    for size in range(1, len(blocks)):
        for chosen, covered in subsets(rows, size):
            if sum(bool(block & chosen) or (block & covered) == block
                   for block in blocks) > size:
                return False
    return True


def verify_weights(rows, denominator, weighted_blocks):
    """Check the weighted bound with integer numerators."""
    weighted_blocks = tuple(weighted_blocks)
    full = (1 << len(rows)) - 1
    if type(denominator) is not int or denominator < 1:
        raise ValueError("denominator must be a positive integer")
    if any(type(block) is not int or not 0 < block <= full
           or type(weight) is not int or weight < 1
           for block, weight in weighted_blocks):
        raise ValueError("expected nonempty blocks and positive integer weights")
    total = sum(weight for _, weight in weighted_blocks)
    for size in range(1, min(len(rows), (total - 1) // denominator) + 1):
        for selected, covered in subsets(rows, size):
            cost = sum(weight for block, weight in weighted_blocks
                       if (block & selected) or (block & covered) == block)
            if cost > size * denominator:
                raise ValueError(f"capacity violated by vertex mask {selected}")
    return total, denominator


def find_reception(rows, count):
    """Enumerate disjoint nonempty blocks with increasing least vertices."""
    tests = [
        (chosen, covered, size)
        for size in range(1, count)
        for chosen, covered in subsets(rows, size)
    ]
    active = [0] * len(tests)
    block_tests = [
        [index for index, (chosen, covered, _) in enumerate(tests)
         if block & chosen or block & covered == block]
        for block in range(1 << len(rows))
    ]

    def search(available, left):
        if not left:
            return ()
        if available.bit_count() < left:
            return None
        first = available & -available
        rest = available ^ first
        tail = rest
        while True:
            block = first | tail
            triggered = block_tests[block]
            if all(active[index] < tests[index][2] for index in triggered):
                for index in triggered:
                    active[index] += 1
                remainder = search(rest ^ tail, left - 1)
                if remainder is not None:
                    return (block,) + remainder
                for index in triggered:
                    active[index] -= 1
            if not tail:
                break
            tail = (tail - 1) & rest
        return search(rest, left)

    return search((1 << len(rows)) - 1, count)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("graphs", nargs="+")
    args = parser.parse_args()
    for record in args.graphs:
        rows = closed_neighborhoods(record)
        full = (1 << len(rows)) - 1
        gamma = next(size for size in range(1, len(rows) + 1)
                     if any(covered == full for _, covered in subsets(rows, size)))
        result = find_reception(rows, gamma)
        if result is not None and not verify_reception(rows, result):
            raise RuntimeError("invalid reception certificate")
        print(record, gamma, result)
