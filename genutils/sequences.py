'''Tools for generating and manipulating ordered sequences of numbers'''

from typing import Generator, Iterable, Sequence, TypeVar, Union

from operator import mul
from collections import defaultdict, Counter

from functools import reduce
from itertools import count, product as cartesian_product


T = TypeVar('T') # generic type for sequence element
U = TypeVar('U') # generic type for a distinct sequence element


def is_unique(seq : Sequence) -> bool:
    '''Whether or not a Sequence contains repeating items'''
    return len(set(seq)) == len(seq)

def product(seq : Iterable[T]) -> T:
    '''Multiplicative analogue to builtin sum()'''
    return reduce(mul, seq)

def int_complement(seq : Sequence[int], bounded : bool=False) -> Generator[int, None, None]:
    '''Generate ordered non-negative integers which don't appear in a collection of integers'''
    _max = max(seq) # cache maximum
    for i in range(_max):
        if i not in seq:
            yield i

    if not bounded: # keep counting past max if unbounded
        yield from count(start=_max + 1, step=1)

def pad_sequence(target_list : Sequence[T], to_length : int, pad_value : U=0, from_left : bool=False) -> list[Union[T, U]]:
    '''Pad a given list with a particular value'''
    padding_list = [pad_value] * (to_length - len(target_list)) # will be empty if the target length is shorter than the provided list (i.e. no padding)
    if from_left:
        return padding_list + target_list
    return target_list + padding_list

def cycle_items(seq : Sequence[T], places : int=1) -> list[T]:
    '''
    Cyclically shift all items in a sequence over by the target number of indices
    By default shifts right, but can also move to left with negative "places" argument
    
    Examples:
        cycle_items([1,2,3,4],  2) -> [3,4,1,2]
        cycle_items([1,2,3,4], -1) -> [4,1,2,3]
    '''
    n_items = len(seq) # this length call is what requires the input to be a Seuquence and not just an Iterable
    return [
        seq[i % n_items]
            for i in range(places, places + n_items)
    ]

def bin_ids_forming_sequence(sequence : Iterable[T], choice_bins : Iterable[Iterable[T]], unique_bins : bool=False) -> Generator[tuple[int, ...], None, None]:
    '''Takes an ordered sequence of N objects and a collection of any number of bins containing arbitrary objects and generates
    all possible N-tuples of bin indices which could produce the target sequence when drawing from those bins WITH replacement
    
    If unique_bins is True, will not return a sequence in which a bin can appear more than once''' 
    occurrences = defaultdict(Counter) # keys are objects, values give indices of bins in which the objects occur, along with multiplicities 
    for i, tags in enumerate(choice_bins):
        for tag in tags:
            occurrences[tag][i] += 1 # NOTE : implementation here requires that T be a hashable type

    for idxs in cartesian_product(*(occurrences[item].keys() for item in sequence)): # TODO : enable no-replacement scheme 
        if (not unique_bins) or is_unique(idxs):
            yield idxs