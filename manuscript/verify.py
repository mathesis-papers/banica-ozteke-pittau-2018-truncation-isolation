import itertools
import json
import math


def is_prime(number):
    return number >= 2 and all(number % divisor for divisor in range(2, math.isqrt(number) + 1))


def field_root(prime):
    multiplier = 1
    while not is_prime(multiplier * prime + 1):
        multiplier += 1
    modulus = multiplier * prime + 1
    for candidate in range(2, modulus):
        root = pow(candidate, multiplier, modulus)
        if root != 1 and pow(root, prime, modulus) == 1:
            return modulus, root
    raise AssertionError('No root of the required order')


def modular_rank(matrix, modulus):
    if not matrix:
        return 0
    rows = [row[:] for row in matrix]
    pivot_index = 0
    for column in range(len(rows[0])):
        selected = next((index for index in range(pivot_index, len(rows)) if rows[index][column] % modulus), None)
        if selected is None:
            continue
        rows[pivot_index], rows[selected] = rows[selected], rows[pivot_index]
        inverse = pow(rows[pivot_index][column], -1, modulus)
        pivot = [(entry * inverse) % modulus for entry in rows[pivot_index]]
        rows[pivot_index] = pivot
        for index in range(pivot_index + 1, len(rows)):
            coefficient = rows[index][column] % modulus
            if coefficient:
                rows[index] = [(entry - coefficient * pivot_entry) % modulus for entry, pivot_entry in zip(rows[index], pivot)]
        pivot_index += 1
        if pivot_index == len(rows):
            break
    return pivot_index


def jacobian_matrix(prime, subset, modulus, root):
    matrix = []
    positions = {vertex: index for index, vertex in enumerate(subset)}
    for first, second in itertools.permutations(subset, 2):
        row = [0] * (len(subset) * prime)
        for column in range(prime):
            coefficient = pow(root, ((first - second) * column) % prime, modulus)
            row[positions[first] * prime + column] = coefficient
            row[positions[second] * prime + column] = -coefficient % modulus
        matrix.append(row)
    return matrix


def check_forest(prime, subset, step):
    selected = set(subset)
    successors = {vertex: (vertex + step) % prime for vertex in subset if (vertex + step) % prime in selected}
    for start in subset:
        current = start
        visited = set()
        while current in successors:
            assert current not in visited
            visited.add(current)
            current = successors[current]
    return len(successors)


def main():
    forest_count = 0
    for prime in [2, 3, 5, 7, 11, 13]:
        for mask in range(1, (1 << prime) - 1):
            subset = [vertex for vertex in range(prime) if mask & (1 << vertex)]
            for step in range(1, prime):
                check_forest(prime, subset, step)
                forest_count += 1
    cases = []
    for prime in [2, 3, 5, 7, 11, 13]:
        modulus, root = field_root(prime)
        subsets = []
        if prime <= 7:
            subsets.extend(tuple(vertex for vertex in range(prime) if mask & (1 << vertex)) for mask in range(1, 1 << prime))
        else:
            subsets.extend(tuple(range(size)) for size in range(1, prime + 1))
            subsets.extend(tuple(vertex for vertex in range(prime) if vertex not in removed) for removed in [(0, 2), (1, 4)])
        for subset in sorted(set(subsets), key=lambda value: (len(value), value)):
            size = len(subset)
            rank = modular_rank(jacobian_matrix(prime, subset, modulus, root), modulus)
            expected = size * (size - 1) if size < prime else (prime - 1) ** 2
            assert rank == expected, (prime, subset, rank, expected)
            defect = size * prime - rank
            dephased_dimension = defect - (size + prime - 1)
            assert dephased_dimension == ((size - 1) * (prime - size - 1) if size < prime else 0)
            cases.append({'p': prime, 'S': list(subset), 'field_modulus': modulus, 'root': root, 'rank': rank, 'defect': defect, 'dephased_dimension': dephased_dimension})
    print(json.dumps({'status': 'pass', 'arithmetic': 'integer modular arithmetic; no floating-point operations', 'scope': 'finite-field Jacobian ranks and induced-cycle forest identities for the listed finite cases', 'forest_checks': forest_count, 'forest_primes': [2, 3, 5, 7, 11, 13], 'jacobian_cases': len(cases), 'cases': cases}, indent=2))


if __name__ == '__main__':
    main()
