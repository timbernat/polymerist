'''Enumeration algorithms for choosing a sequence of symbols out of a sequence of sets of those symbols'''

from typing import Generator, Iterable, Sequence, TypeVar
T = TypeVar('T') # generic type for sequence element

from copy import deepcopy
from collections import defaultdict, Counter
from itertools import product as cartesian_product

from .seqops import is_unique


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