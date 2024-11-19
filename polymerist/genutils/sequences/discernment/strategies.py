'''Abstract base and concrete implementations of algorithms which solve the DISCERNMENT Problem'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Generator, Generic, ParamSpec, Sequence, TypeVar
from abc import ABC, abstractmethod

from copy import deepcopy # TODO : deprecate the need for this in Cartesian
from itertools import product as cartesian_product

P = ParamSpec('P')
from .inventory import SymbolInventory, T, L


# HELPER FUNCTIONS
def is_unique(seq : Sequence) -> bool:
    '''Whether or not a Sequence contains repeating items'''
    return len(set(seq)) == len(seq)

# ABSTRACT BASE WHICH DEFINES STRATEGY INTERFACE
class DISCERNMENTStrategy(ABC, Generic[T, L]):
    '''
    DISCERNMENT = the Determination of Index Sequences from Complete Enumeration of Ransom Notes (Multiset Extension with Nonlexical Types)
    Provides an implementation of the general ransom note enumeration problem
    '''
    @abstractmethod
    def enumerate_choice_labels(self,
        word : Sequence[T],
        symbol_inventory : SymbolInventory[T, L],
        ignore_multiplicities : bool=False,
        unique_bins : bool=False,
        **kwargs : P.kwargs,
    ) -> Generator[tuple[L, ...], None, None]:
        '''
        Takes a word (N-element sequence of type T) and a symbol inventory (map from symbols and bin labels to counts),
        Exhaustively enumerates all N-tuples of indices corresponding to ordering of bins from which the word could be drawn
        
        Support modifications to base behavior:
        *** If ignore_multiplicities=True, will not respect the counts of elements in each bin when drawing
            For a given symbol, this allows a bin containing the symbol to appear more times than that symbol is present in that bin
        *** If unique_bins=True, will only allow each bin to be sampled from once, 
            EVEN if that bin contains symbols appearing later in the word
        '''
        pass

# CONCRETE IMPLEMENTATIONS
class DISCERNMENTStrategyCartesian(DISCERNMENTStrategy):
    '''
    Naive implementation of generalized ransom note enumeration strategy based on cartesian products
    
    Generates all possible index sequences ignoring multiplicity, then checks them one-by-one to see if they're valid

    Roughly 1-2 OOM slower (prefactor) than Recursive and Stack
    unless ignoring multiplicity in which case this is ~1 OOM faster
    '''
    def enumerate_choice_labels(self,
        word : Sequence[T],
        symbol_inventory : SymbolInventory[T, L],
        ignore_multiplicities : bool=False,
        unique_bins : bool=False,
    ) -> Generator[tuple[L, ...], None, None]:
        # generate every valid bin index tuple WITHOUT regard to repetition or uniqueness
        for indices in cartesian_product(*(symbol_inventory[symbol].keys() for symbol in word)): 
            if unique_bins and not is_unique(indices): 
                continue # skip non-unique bin choices if the option is set

            if not ignore_multiplicities:
                choice_inventory = deepcopy(symbol_inventory) # make a new copy for each path check
                for i, sym in zip(indices, word, strict=True):
                    if choice_inventory[sym][i] == 0:
                        overdrawn = True
                        break
                    choice_inventory[sym][i] -= 1
                else:
                    overdrawn = False 

                if overdrawn: # TOSELF: there are some optimizations to be made here (i.e. avoid copying inventory via backtrack)...
                    continue  # ... but am keeping this exact implementation for backward-compatibility reasons
            
            yield indices # only yield if all specified conditions are met

class DISCERNMENTStrategyRecursive(DISCERNMENTStrategy):
    '''
    Recursive implementation of generalized ransom note enumeration strategy via divide-and-conquer approach

    Breaks sequence down into first and all remaining symbols ("head" and "tail", respectively)
    Yields solution as all indices where the head occurs + recursive solutions to the tail sequence enumeration

    Easiest to analyze and reason about, but imposes cap on word length due to Python call stack size
    Performs intermediately in benchmark (faster than Cartesian, but slower than Stack)
    '''
    def enumerate_choice_labels(self,
        word : Sequence[T],
        symbol_inventory : SymbolInventory[T, L],
        ignore_multiplicities : bool=False,
        unique_bins : bool=False,
        _buffer : tuple[int]=None,
    ) -> Generator[tuple[L, ...], None, None]:
        if _buffer is None:
            _buffer = tuple()

        if not word: # base case for recursion
            yield _buffer # yields empty buffer if no sequence is present
            return

        symbol_cost = (0 if ignore_multiplicities else 1) # NOTE: definition here could be made more terse, but at the expense of readability
        
        symbol, *tail = word # separate head symbol from rest of sequence
        for bin_idx, num_occurences in symbol_inventory[symbol].items():
            # only proceed if letters are available from that bin, AND either uniqueness is not required, or it is required but has not yet been violated
            if (num_occurences > 0) and (not unique_bins or (bin_idx not in _buffer)):
                symbol_inventory[symbol][bin_idx] -= symbol_cost # mark current symbol and bin in symbol inventory and visit tracker, respectively

                yield from self.enumerate_choice_labels( # recursive traversal of remainder of sequence with current choice of bin for leading character
                    word=tail,
                    symbol_inventory=symbol_inventory,
                    ignore_multiplicities=ignore_multiplicities,
                    unique_bins=unique_bins,
                    _buffer=_buffer + (bin_idx,), # creates copy, rather than modifying the buffer for the current symbol (i.e. exactly what we want)
                )
                symbol_inventory[symbol][bin_idx] += symbol_cost # replace removed symbol for subsequent traversals to avoid polluting other states

class DISCERNMENTStrategyStack(DISCERNMENTStrategy):
    '''
    Stack-based implementation of generalized ransom note enumeration strategy

    Pushes all valid symbol-index pairs onto a stack for the current word positions,
    then advances and does the same until reaching the final symbol, backtracking and restoring symbols afterwards

    Fastest of all strategies across benchmark (except Cartesian **specifically** when ignore_multiplicities=True)
    '''
    def enumerate_choice_labels(self,
        word : Sequence[T],
        symbol_inventory : SymbolInventory[T, L],
        ignore_multiplicities : bool=False,
        unique_bins : bool=False,
    ) -> Generator[tuple[L, ...], None, None]:
        symbol_cost = (0 if ignore_multiplicities else 1) # NOTE: definition here could be made more terse, but at the expense of readability

        buffer = []
        bin_stack : list[tuple[T, int, bool]] = [(None, None, False)] # initialize with null symbol and bin, with no intent to replace (i.e. with intent to take)

        while bin_stack:
            symbol, bin_id, replace = bin_stack.pop()
            if replace:
                _last_bin_id = buffer.pop() # NOTE: this happens to be the same as bin_id (by construction) but will keep separate for clarity
                symbol_inventory[symbol][bin_id] += symbol_cost
            else:   
                if bin_id is not None: # exclude initial stack to give correct behavior in the case of an empty sequence (no explicit base case required!)
                    symbol_inventory[symbol][bin_id] -= symbol_cost
                    buffer.append(bin_id)

                if len(buffer) == len(word): # report buffer if the whole sequence has been accounted for
                    yield tuple(buffer)
                else: 
                    next_symbol = word[len(buffer)] # NOTE: symbol update fetches need to be done carefully WRT buffer updates
                    for next_bin_id, count in reversed(symbol_inventory[next_symbol].items()): # NOTE : reversal is solely to ensure overall order is consistent with other methods (needed since a stack is a FILO queue)
                        if (count > 0) and (not unique_bins or (next_bin_id not in buffer)):
                            bin_stack.append( (next_symbol, next_bin_id, True) )  # push replace operation onto stack first so it happens last
                            bin_stack.append( (next_symbol, next_bin_id, False) ) # presence of single take and replace