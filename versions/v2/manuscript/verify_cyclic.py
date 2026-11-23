import itertools
import json
import math

from verify import is_prime, jacobian_matrix, modular_rank


def prime_factors(number):
    return [factor for factor in range(2, number + 1) if number % factor == 0 and is_prime(factor)]


def exact_root(order):
    multiplier = 1
    while not is_prime(multiplier * order + 1):
        multiplier += 1
    field_modulus = multiplier * order + 1
    for candidate in range(2, field_modulus):
        root = pow(candidate, multiplier, field_modulus)
        if all(pow(root, order // factor, field_modulus) != 1 for factor in prime_factors(order)):
            return field_modulus, root
    raise AssertionError('No root of the required exact order')


def cycle_count(order, subset):
    selected = set(subset)
    total = 0
    for step in range(1, order):
        divisor = math.gcd(order, step)
        total += sum(all(vertex in selected for vertex in range(residue, order, divisor))
                     for residue in range(divisor))
    return total


def component_defect(order, subset):
    total = len(subset)
    for step in range(1, order):
        remaining = set(subset)
        components = 0
        while remaining:
            pending = [remaining.pop()]
            components += 1
            while pending:
                vertex = pending.pop()
                for neighbor in [(vertex + step) % order, (vertex - step) % order]:
                    if neighbor in remaining:
                        remaining.remove(neighbor)
                        pending.append(neighbor)
        total += components
    return total


def abelian_check(dimensions, subset):
    group = list(itertools.product(*(range(dimension) for dimension in dimensions)))
    order = len(group)
    exponent = math.lcm(*dimensions)
    field_modulus, root = exact_root(exponent)
    positions = {vertex: index for index, vertex in enumerate(subset)}
    matrix = []
    for first, second in itertools.permutations(subset, 2):
        row = [0] * (len(subset) * order)
        for column_index, character in enumerate(group):
            power = sum((first[index] - second[index]) * character[index] * (exponent // dimension)
                        for index, dimension in enumerate(dimensions))
            coefficient = pow(root, power % exponent, field_modulus)
            row[positions[first] * order + column_index] = coefficient
            row[positions[second] * order + column_index] = -coefficient % field_modulus
        matrix.append(row)
    rank = modular_rank(matrix, field_modulus)
    complete_cosets = 0
    component_total = len(subset)
    selected = set(subset)
    for step in group:
        if not any(step):
            continue
        remaining = set(group)
        while remaining:
            seed = next(iter(remaining))
            orbit = set()
            vertex = seed
            while vertex not in orbit:
                orbit.add(vertex)
                vertex = tuple((entry + shift) % dimension
                               for entry, shift, dimension in zip(vertex, step, dimensions))
            remaining.difference_update(orbit)
            complete_cosets += orbit <= selected
        remaining = set(subset)
        while remaining:
            pending = [remaining.pop()]
            component_total += 1
            while pending:
                vertex = pending.pop()
                for sign in [-1, 1]:
                    neighbor = tuple((entry + sign * shift) % dimension
                                     for entry, shift, dimension in zip(vertex, step, dimensions))
                    if neighbor in remaining:
                        remaining.remove(neighbor)
                        pending.append(neighbor)
    expected_defect = len(subset) * (order - len(subset) + 1) + complete_cosets
    assert len(subset) * order - rank == component_total == expected_defect
    return {'cyclic_factors': list(dimensions), 'S': [list(vertex) for vertex in subset],
            'field_modulus': field_modulus, 'root_order': exponent, 'root': root,
            'complete_step_cosets': complete_cosets, 'rank': rank, 'defect': expected_defect}


def main():
    combinatorial_cases = 0
    for order in range(2, 13):
        for mask in range(1, 1 << order):
            subset = tuple(vertex for vertex in range(order) if mask & (1 << vertex))
            expected = len(subset) * (order - len(subset) + 1) + cycle_count(order, subset)
            assert component_defect(order, subset) == expected
            combinatorial_cases += 1
    cases = []
    for order in [2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15]:
        if order <= 8:
            subsets = {tuple(vertex for vertex in range(order) if mask & (1 << vertex))
                       for mask in range(1, 1 << order)}
        else:
            subsets = {tuple(range(size)) for size in range(1, order + 1)}
            subsets.update(tuple(range(0, order, divisor))
                           for divisor in range(2, order) if order % divisor == 0)
            subsets.add(tuple(vertex for vertex in range(order) if vertex not in [0, 2]))
        field_modulus, root = exact_root(order)
        for subset in sorted(subsets, key=lambda value: (len(value), value)):
            size = len(subset)
            cycles = cycle_count(order, subset)
            rank = modular_rank(jacobian_matrix(order, subset, field_modulus, root), field_modulus)
            expected_rank = size * (size - 1) - cycles
            assert rank == expected_rank, (order, subset, rank, expected_rank)
            defect = size * order - rank
            assert defect == component_defect(order, subset)
            dephased_defect = defect - (size + order - 1)
            assert dephased_defect >= 0
            cases.append({'N': order, 'S': list(subset), 'field_modulus': field_modulus,
                          'root': root, 'complete_step_cosets': cycles, 'rank': rank,
                          'defect': defect, 'dephased_defect': dephased_defect,
                          'orthogonality_map_regular': cycles == 0})
    abelian_cases = []
    for dimensions in [(2, 2), (2, 3), (2, 4), (2, 2, 2), (3, 3)]:
        group = list(itertools.product(*(range(dimension) for dimension in dimensions)))
        order = len(group)
        masks = range(1, 1 << order) if order <= 8 else [(1 << size) - 1 for size in range(1, order + 1)]
        for mask in masks:
            subset = tuple(vertex for index, vertex in enumerate(group) if mask & (1 << index))
            abelian_cases.append(abelian_check(dimensions, subset))
    print(json.dumps({'status': 'pass', 'arithmetic': 'integer modular arithmetic; no floating-point operations',
                      'scope': 'finite abelian Fourier tangent dimensions; finite ranks do not assert integration at singular points',
                      'combinatorial_cases': combinatorial_cases,
                      'all_subsets_combinatorial_orders': list(range(2, 13)),
                      'jacobian_cases': len(cases), 'cases': cases,
                      'abelian_jacobian_cases': len(abelian_cases), 'abelian_cases': abelian_cases}, indent=2))


if __name__ == '__main__':
    main()
