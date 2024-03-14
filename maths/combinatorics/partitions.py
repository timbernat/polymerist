'''For explicitly enumerating partitions of sets and multisets'''

from typing import Generator

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