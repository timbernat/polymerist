'''Utilites for representing pure permutations, cycles, and permutation groups'''

from typing import Generator, Iterable, Optional, Sequence, TypeVar
from dataclasses import dataclass, field
T = TypeVar('T')

import numpy as np
from math import lcm
from operator import mul

from functools import reduce
from itertools import chain, permutations
from collections import Counter, defaultdict


def int_to_factoradic(n : int) -> list[int]:
    '''Determine the digits of the factorial base representation of an integer'''
    if n < 0:
        raise ValueError('Cannot generate factoradic digits of a negative integer')
    if n == 0:
        return [0]

    digits : list[int] = []
    quotient, divisor = n, 1
    while quotient > 0:
        quotient, digit = divmod(quotient, divisor)
        digits.append(digit)
        divisor += 1
    return digits[::-1]

def factoradic_to_int(digits : Sequence[int]) -> int:
    '''Convert the digits of the factorial base representation of an integer back to a standard base-10 integer'''
    n = 0
    place = 1
    for i, digit in enumerate(reversed(digits), start=1):
        if not (0 <= digit <= i):
            raise ValueError(f'Digit {digit} in {place}!-place exceeds max allowed factoradic value of {i} for that place')
        n += place * digit
        place *= i
    return n


class Cycle(tuple):
    '''For representing a cyclic collection of objects'''
    def __new__(cls, *values : Iterable[T]) -> None: # override new to allow for star unpacking
        return super(Cycle, cls).__new__(cls, values)

    def __repr__(self) -> str:
        '''Represent as a '''
        return f'{self.__class__.__name__}{super().__repr__()}'

    def __getitem__(self, __key : int) -> T:
        '''Access an item, with wrap-around indices'''
        return super().__getitem__(__key % len(self))

    def __reversed__(self) -> 'Cycle':
        return self.__class__(*super().__reversed__())

    def starting_from_index(self, index : int) -> 'Cycle':
        '''Return a Cycle which is in the same order as self but begins at the given index'''
        return Cycle(*(self[i] for i in range(index, index + len(self))))

    def __add__(self, other : 'Cycle') -> 'Cycle':
        '''Extend a Cycle by another Cycle'''
        if not isinstance(other, Cycle):
            raise TypeError
        return self.__class__(*super().__add__(other))

    def __call__(self, other) -> 'Cycle': # TODO : implement composition via R-to-L adjacency
        ...

    def copy(self) -> 'Cycle':
        return self.__class__(*self)

    @property
    def mapping(self) -> dict[T, T]:
        '''A mapping of each item in the cycle to the following item'''
        return {
            self[i - 1] : elem # NOTE: implemented this way specifically to make use of "-1" last index and avoid and IndexError for the wrap-around
                for i, elem in enumerate(self)
        }

    # HELPER METHODS
    @staticmethod
    def max_elem_in_cycles(cycles : Iterable['Cycle']) -> int:
        '''Returns the largest element across a collection of Cycles'''
        return max(chain(*cycles))

    @staticmethod
    def cycles_are_disjoint(cycles : Iterable['Cycle']) -> bool:
        '''Check if a cycle decomposition is disjoint'''
        return set.intersection(*map(set, cycles)) == set() # check that the intersection of all cycles is the empty set

    @staticmethod
    def cycles_produce_partition(cycles : Iterable['Cycle']) -> bool:
        '''Check if a cycle decomposition forms a partition of some set of integers'''
        degree = Cycle.max_elem_in_cycles(cycles) + 1 # find total maximum to use as order
        all_elems = set.union(*map(set, cycles)) # check that the intersection of all cycles is the empty set

        for i in range(degree):
            if i not in all_elems:
                return False
        else:
            return True

    @staticmethod
    def cycle_type(cycles : Iterable['Cycle']) -> dict[int, int]:
        '''Returns a dict of cycle lengths and the number of cycle in the permutation with that length'''
        cycle_len_counts = Counter()
        for cycle in cycles:
            cycle_len_counts[len(cycle)] += 1
        longest_cycle_len = max(cycle_len_counts.keys())

        return {
            cycle_len : cycle_len_counts[cycle_len]
                for cycle_len in range(1, longest_cycle_len + 1) # include zeros for intermediate sizes, 1-index to make use of counts
        }

    @staticmethod
    def cycle_index_sym(cycle_type : dict[int, int], var_sym : str='x', mul_sym : str='*', exp_sym : str='^') -> str:
        '''Represent a single cycle type in symbolic variable form
        Can optionally supply custom symbols for variable (single-char only), multiplication sign, and exponentiation sign'''
        if len(var_sym) > 1:
            raise ValueError('must provide sinle-character symbolic variable for cycle index')
        return mul_sym.join(f'{var_sym}_{cycle_len}{exp_sym}{count}' for cycle_len, count in cycle_type.items())
