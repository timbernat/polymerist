'''Tools for generating and manipulating ordered sequences of numbers'''

from typing import Generator, Iterable, TypeVar

from operator import mul
from functools import reduce
from itertools import count


T = TypeVar('T') # generic type for sequence element
def product(seq : Iterable[T]) -> T:
    '''Analogous to builtin sum()'''
    return reduce(mul, seq)

def int_complement(seq : Iterable[int], bounded : bool=False) -> Generator[int, None, None]:
    '''Generate ordered non-negative integers which don't appear in a collection of integers'''
    _max = max(seq) # cache maximum
    for i in range(_max):
        if i not in seq:
            yield i

    if not bounded: # keep counting past max if unbounded
        yield from count(start=_max + 1, step=1)