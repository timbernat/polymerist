'''Tools for generating and manipulating ordered sequences of numbers'''

from typing import Generator, Iterable, TypeVar

from operator import mul
from collections import defaultdict, Counter

from functools import reduce
from itertools import count, product as cartesian_product


T = TypeVar('T') # generic type for sequence element
def product(seq : Iterable[T]) -> T:
    '''Multiplicative analogue to builtin sum()'''
    return reduce(mul, seq)

def int_complement(seq : Iterable[int], bounded : bool=False) -> Generator[int, None, None]:
    '''Generate ordered non-negative integers which don't appear in a collection of integers'''
    _max = max(seq) # cache maximum
    for i in range(_max):
        if i not in seq:
            yield i

    if not bounded: # keep counting past max if unbounded
        yield from count(start=_max + 1, step=1)

def bin_ids_forming_sequence(sequence : Iterable[T], choice_bins : Iterable[Iterable[T]]) -> Generator[tuple[int, ...], None, None]:
    '''Takes an ordered sequence of N objects and a collection of any number of bins containing arbitrary objects and generates
    all possible N-tuples of bin indices which could produce the target sequence when drawing from those bins WITH replacement'''
    occurrences = defaultdict(Counter) # keys are objects, values give indices of bins in which the objects occur, along with multiplicities 
    for i, tags in enumerate(choice_bins):
        for tag in tags:
            occurrences[tag][i] += 1 # NOTE : implementation here requires that T be a hashable type

    for idxs in cartesian_product(*(occurrences[item].keys() for item in sequence)): # TODO : enable no-replacement scheme 
        yield idxs