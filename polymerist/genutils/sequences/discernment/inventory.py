'''Utilities and type-hinting for creating symbol inventories, which map symbols and bin labels to occurences
Useful concrete data structure for representing ordered sequences of symbol multisets for generalized ransom note enumeration'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import logging
LOGGER = logging.getLogger(__name__)

from typing import Any, Generator, Generic, Iterable, Mapping, Optional, Sequence, TypeVar, Union

T = TypeVar('T')
L = TypeVar('L')

from importlib.util import find_spec
from functools import cached_property
from collections import Counter, defaultdict

try:
    import numpy as np
except ModuleNotFoundError:
    pass


def full_arr_builtin(*dims : Iterable[int], fill_value : Any=None) -> list[Any]:
    '''Similar to numpy.full (https://numpy.org/doc/stable/reference/generated/numpy.full.html) but with builtin list type'''
    arrlist = fill_value
    for dim_size in reversed(dims):
        if not isinstance(dim_size, int):
            raise TypeError(f'All dimensions must be of type int (not {dim_size.__class__.__name__})')
        arrlist = [arrlist for _ in range(dim_size)]
    return arrlist

class SymbolInventory(dict, Generic[T, L]):
    '''
    Representation of a map from a set of symbols of type T and a set of bin labels of type L to integer counts
    Implemented as a dict, keyed by symbol, whose values count the number of occurrences of that symbol in all bins that symbol occurs

    Data structure with specific methods which are useful in solving the generalized ransom note enumeration problem
    '''
    # INSTANTIATION
    def __init__(
        self, *args,
         _number_of_symbols : Optional[int]=None,
        _number_of_bins     : Optional[int]=None,
        _symbol_index_map   : Optional[dict[T, int]]=None,
        _bin_index_map      : Optional[dict[L, int]]=None,
        **kwargs,
    ) -> None:
        self._number_of_symbols = _number_of_symbols
        self._number_of_bins    = _number_of_bins

        self._symbol_index_map  = _symbol_index_map
        self._bin_index_map     = _bin_index_map

        super().__init__(*args, **kwargs)

    @classmethod
    def from_bin_sequence_mapping(cls,
        choice_bin_map : Mapping[L, Sequence[Iterable[T]]],
        _symbol_index_map : Optional[dict[T, int]]=None,
        _bin_index_map    : Optional[dict[L, int]]=None,
    ) -> 'SymbolInventory[T, L]':
        '''Initialize inventory from a mapping of labels to bins'''
        sym_inv_dict = defaultdict(Counter[L, int]) # keys are objects of type T ("symbols"), values give multiplicities of symbols keyed by bin position
        for label, choice_bin in choice_bin_map.items():
            for sym in choice_bin:
                sym_inv_dict[sym][label] += 1 # NOTE : implementation here requires that T be a hashable type
        
        return cls(
            sym_inv_dict,
            _number_of_symbols=len(sym_inv_dict), # pre-compute numbers of symbols and bins
            _number_of_bins =len(choice_bin_map), # pre-compute numbers of symbols and bins
            _symbol_index_map=_symbol_index_map, # enable iltering down of pre-computed maps from other initialization methods
            _bin_index_map=_bin_index_map,       # enable iltering down of pre-computed maps from other initialization methods
        )
    
    @classmethod
    def from_bin_sequence(cls, choice_bins : Sequence[Iterable[T]]) -> 'SymbolInventory[T, int]':
        '''Initialize inventory from an ordered sequence of bins
        Special case of mapping instantiation for sequences of bins; simply uses the index of a bin as its label'''

        choice_bin_map : dict[int, Sequence[Iterable[T]]] = {}
        _bin_index_map : dict[int, int] = {}
        for i, choice_bin in enumerate(choice_bins):
            choice_bin_map[i] = choice_bin # expand sequence into dict keyed by position in sequence 
            _bin_index_map[i] = i # ensures bin order does not get scrambled when initializing from this route (avoids confusion down the line)

        return cls.from_bin_sequence_mapping(
            choice_bin_map=choice_bin_map,
            _bin_index_map = _bin_index_map,
        )
    
    @classmethod
    def from_bins(cls, choice_bins : Union[Sequence[Iterable[T]], dict[L, Sequence[Iterable[T]]]]) -> 'SymbolInventory[T, L]':
        '''"Smart" initialization method which dispatches to from_bin_sequence() or from_bin_sequence_mapping(),
        depending on the nature of the "choice_bins" container provided'''
        if isinstance(choice_bins, Mapping):
            return cls.from_bin_sequence_mapping(choice_bin_map=choice_bins)
        elif isinstance(choice_bins, Sequence):
            return cls.from_bin_sequence(choice_bins=choice_bins)
        elif isinstance(choice_bins, Generator):
            return cls.from_bin_sequence(choice_bins=[i for i in choice_bins])
        else:
            raise TypeError(f'{cls.__name__} does not support initialization from non Sequence/Mapping container of type "{choice_bins.__class__.__name__}"')

    # INSTANCE METHODS
    def __repr__(self) -> str:
        '''Wrap dict repr to clarify that this is a subclass'''
        return f'{self.__class__.__name__}({super().__repr__()})'
    
    @property
    def number_of_symbols(self) -> int:
        '''Number of unique symbols present in the SymbolInventory'''
        if self._number_of_symbols is None:
            LOGGER.debug(f'Initializing field {"number_of_symbols"} of {self!r}')
            self._number_of_symbols = len(self.keys())
        return self._number_of_symbols
    num_symbols = num_sym = n_symbols = n_sym = number_of_symbols # aliases for convenience 

    @property
    def number_of_bins(self) -> int:
        '''Number of unique symbols present in the SymbolInventory'''
        if self._number_of_bins is None:
            LOGGER.debug(f'Initializing field {"number_of_bins"} of {self!r}')
            self._number_of_bins = self.involution.number_of_symbols
        return self._number_of_bins
    num_bins = num_bin = n_bins = n_bin = number_of_bins # aliases for convenience

    @property
    def symbol_index_map(self) -> dict[T, int]:
        '''Arbitrary one-to-one mapping between all symbols present in the inventory and integer indices'''
        if self._symbol_index_map is None:
            LOGGER.debug(f'Initializing field {"symbol_index_map"} of {self!r}')
            self._symbol_index_map = {
                symbol : i
                    for i, symbol in enumerate(self.keys())
            }
        return self._symbol_index_map

    @property
    def bin_index_map(self) -> dict[T, int]:
        '''Arbitrary one-to-one mapping between all bins present in the inventory and integer indices'''
        if self._bin_index_map is None:
            LOGGER.debug(f'Initializing field {"bin_index_map"} of {self!r}')
            self._bin_index_map = self.involution.symbol_index_map
        return self._bin_index_map
    
    def contains_word(self, sequence : Sequence[T], ignore_multiplicities : bool=False) -> bool:
        '''
        Check if a word (i.e. sequence of symbols of type T) could possibly be produced from a SymbolInventory

        Returns bool; True indicated the word can be made from the inventory
        False return indicates the sequence contains symbols not present in the inventory 
        OR contains more of one particular symbol than are present in the whole inventory
        '''
        for symbol, count in Counter(sequence).items():
            if symbol not in self:
                return False
            elif self[symbol].total() < (1 if ignore_multiplicities else count): # at least one occurence suffices if we don't care about exact counts
                return False # NOTE: could merge with above check (via OR operator short-circuit), but separated here for readability
        else:
            return True

    ## CONVERSION TO ALTERNATE FORMS
    def deepcopy(self) -> 'SymbolInventory[T, L]': # NOTE : deliberately named this "deepcopy" to reserve the name for a potetial cheaper shallow-copy "copy" method
        '''Create a deep copy of the current SymbolInventory'''
        return self.__class__(
            {symbol : bin_counts.copy() for symbol, bin_counts in self.items()}, # TOSELF : run some more tests to verify Counter.copy() does not point to original Counter
            _number_of_symbols=self._number_of_symbols,
            _number_of_bins=self._number_of_bins,
            _symbol_index_map=self._symbol_index_map,
            _bin_index_map=self._bin_index_map,
        )

    @cached_property
    def involution(self) -> 'SymbolInventory[L, T]':
        '''Returns a new SymbolInventory which switches the order of precedence of the symbols and bin labels
        So-called because self.involution.involution returns a SymbolInventory identical to self (by construction)'''
        invol = defaultdict(Counter[T, int])
        for symbol, bin_counts in self.items():
            for bin_label, count in bin_counts.items():
                invol[bin_label][symbol] += count
        return self.__class__(invol)
    
    @cached_property
    def occurence_matrix(self) -> Union[list[list[int]], np.ndarray[Any, int]]:
        '''
        Creates an NxM matrix (N is the number of symbols, and M is the number of bins) which represents a SymbolInventory
        Matrix element A_ij denotes the number of occurences of the i-th symbol (according to an arbitrary numbering) in the j-th bin

        Will attempt to return matrix as a 2-D numpy array, but will default to a doubly-nested list if numpy is not found to be installed
        '''
        shape : tuple[int, int] = (self.number_of_symbols, self.number_of_bins)
        if find_spec('numpy') is not None:
            occ_matr = np.zeros(shape, dtype=int)
        else: # default to list if numpy is not installed
            occ_matr = full_arr_builtin(*shape, fill_value=0)
        
        for symbol, bin_counts in self.items():
            i = self.symbol_index_map[symbol]
            for bin_label, count in bin_counts.items():
                j = self.bin_index_map[bin_label]
                occ_matr[i, j] = count
        
        return occ_matr