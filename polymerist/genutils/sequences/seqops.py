'''Generic operations for indexing, generating, and iterating over sequences'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Generator, Sequence, TypeVar, Union
T = TypeVar('T') # generic type for sequence element
S = TypeVar('S') # generic type for a distinct sequence element
from itertools import count


def is_unique(seq : Sequence) -> bool:
    '''Whether or not a Sequence contains repeating items'''
    return len(set(seq)) == len(seq)

def int_complement(integers : Sequence[int], bounded : bool=False) -> Generator[int, None, None]:
    '''Generate ordered non-negative integers which don't appear in a sequence of integers'''
    _max = max(integers) # cache maximum (precludes use of generator-like sequence)
    for i in range(_max):
        if i not in integers:
            yield i

    if not bounded: # keep counting past max if unbounded
        yield from count(start=_max + 1, step=1)

def pad_sequence(target_list : Sequence[T], to_length : int, pad_value : S=0, from_left : bool=False) -> list[Union[T, S]]:
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
    n_items = len(seq) # this length call is what requires the input to be a Sequence and not just an Iterable
    return [
        seq[i % n_items]
            for i in range(places, places + n_items)
    ]