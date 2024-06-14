'''Generic operations for indexing, generating, and iterating over sequences'''

from typing import Generator, Iterable, Sequence, TypeVar, Union
T = TypeVar('T') # generic type for sequence element
S = TypeVar('S') # generic type for a distinct sequence element

from copy import deepcopy
from collections import defaultdict, Counter
from itertools import count, product as cartesian_product


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
    n_items = len(seq) # this length call is what requires the input to be a Seuquence and not just an Iterable
    return [
        seq[i % n_items]
            for i in range(places, places + n_items)
    ]

def bin_ids_forming_sequence(sequence : Sequence[T], choice_bins : Sequence[Iterable[T]], draw_without_repeats : bool=True, unique_bins : bool=False) -> Generator[tuple[int, ...], None, None]:
    '''
    Takes an ordered sequence of N objects of a given type and an ordered of any number of bins, each containing an arbitary amount of unordered objects of the same type
    Generates all possible N-tuples of bin indices which could produce the target sequence when drawing from those bins in the 
    
    if draw_without_repeats=True, will respect the multiplicity of elements in each bin when drawing
    (i.e. will never have a bin position appear for a given object more times that that object appears in the corresponding bin)

    if unique_bins=True, will only allow each bin to be sampled from once, EVEN if that bin contains elements which may occur later in the sequence
    '''
    symbol_inventory = defaultdict(Counter) # keys are objects of type T ("symbols"), values give multiplicities of symbols keyed by bin position
    for i, choice_bin in enumerate(choice_bins):
        for sym in choice_bin:
            symbol_inventory[sym][i] += 1 # NOTE : implementation here requires that T be a hashable type

    for idxs in cartesian_product(*(symbol_inventory[item].keys() for item in sequence)): # generate every valid sequence of bin positions WITHOUT regard to repetition or uniqueness
        if unique_bins and not is_unique(idxs):
            continue # skip non-unique bin choices if the option is set

        if draw_without_repeats:
            choice_inventory = deepcopy(symbol_inventory) # make a new copy for each path check
            for i, sym in zip(idxs, sequence, strict=True):
                if choice_inventory[sym][i] == 0:
                    overdrawn = True
                    break
                choice_inventory[sym][i] -= 1
            else:
                overdrawn = False

            if overdrawn:
                continue
        
        yield idxs # only yield if all specified conditions are met