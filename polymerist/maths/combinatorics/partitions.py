'''For explicitly enumerating partitions of sets and multisets'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Generator, Sequence


def int_partitions(n : int, _pivot : int=1) -> Generator[tuple[int], None, None]:
    '''Enumerates all integer partitions of n'''
    yield (n,)
    for i in range(_pivot, n//2 + 1):
        for p in int_partitions(n - i, _pivot=i):
            yield (i,) + p

def multiset_partition(n : int, k : int) -> Generator[tuple[int], None, None]:
    '''Enumerates all multisets of k integers whose sum is n (respects order and includes 0s)'''
    if k == 0:
        yield ()
    elif k == 1:
        yield (n,)
    else:
        for i in range(0, n + 1):
            for subpart in multiset_partition(n - i , k - 1):
                yield (i,) + subpart

def make_change_greedy(n : int, denoms : Sequence[int]) -> dict[int, int]:
    '''Greedy algorithm for making change for a given total'''
    # precheck for uniqueness and integrality of denominations
    seen = set()
    for den in denoms:
        if den in seen:
            raise ValueError('No duplicate denominations allowed')
        elif not isinstance(den, int):
            raise TypeError('All denominations must be integers')
        else:
            seen.add(den)

    # actual change-making happens here
    multips = {}
    for den in sorted(denoms, reverse=True):
        d, n = divmod(n, den)
        multips[den] = d

    if n != 0:
        raise ValueError(f'Could not successfully make complete change ({n} left over)')
    return multips