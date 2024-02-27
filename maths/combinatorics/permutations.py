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


@dataclass # this is purely to give nice bool, repr, and eq boilerplates
class Permutation:
    '''For representing a permutation'''
    elems : Sequence[int] = field(default_factory=list)

    def __init__(self, *elems : Sequence[int]) -> None: # need to implement __init__ in order to allow star unpacking
        self.elems = elems
        if self.element_set != set(self.natural_order):
            raise ValueError

    # CARDINALITY
    @property
    def degree(self) -> int:
        return len(self.elems)

    def __len__(self) -> int:
        return self.degree
    
    @property
    def element_set(self) -> set[int]:
        return set(self.elems)
    
    @staticmethod
    def _natural_order(N : int) -> list[int]:
        '''The first N natural numbers {0, ... ,N-1} in ascending order'''
        return [i for i in range(N)]

    @classmethod
    def from_degree(cls, degree : int, random : bool=False) -> 'Permutation':
        '''Construct a permutation of the given degree. Is just the identity by default, but can be made random with the "random" flag'''
        elems = cls._natural_order(degree)
        if random:
            np.random.shuffle(elems)
        return cls(*elems)
    
    @classmethod
    def identity(cls, degree : int) -> 'Permutation':
        '''Construct the Identity permutation of a given degree'''
        return cls.from_degree(degree, random=False)

    # ELEMENTARY PERMUTATIONS
    @property
    def natural_order(self) -> list[int]:
        return self._natural_order(self.degree)
    
    @property
    def as_identity(self) -> 'Permutation':
        '''The identity permutation with the same degree as the current permutation'''
        return self.__class__.identity(self.degree)

    @property
    def inverse(self) -> 'Permutation':
        '''The permutation which, when composed with this permutation from either the left or the right, yields the identity permutation'''
        return self.__class__(*np.argsort(self.elems))
    
    @property
    def reverse(self) -> 'Permutation':
        '''The reversed-order of a permutation'''
        return self.__class__(*reversed(self.elems))
    
    # COMPOSITIONS AND MAPS
    def image(self, coll : Iterable[T]) -> Iterable[T]:
        '''The image of an ordered collection under the defined permutation'''
        # if len(coll) != self.degree:
        #     raise ValueError
        return [coll[i] for i in self.elems]
    
    def __call__(self, elem : int) -> int:
        '''Return the image of a single element under the permutation'''
        if not isinstance(elem, int):
            raise TypeError
        
        if (elem < 0) or (elem > (self.degree - 1)):
            raise ValueError(f'Integer {elem} has no image under permutation of degree {self.degree}')
        
        return self.elems[elem]
    
    def __getitem__(self, elem : int) -> int:
        return self.elems[elem]

    def __mul__(self, other : 'Permutation') -> 'Permutation':
        if not isinstance(other, Permutation):
            raise TypeError(f'Cannot compose {self.__class__.__name__} with {other.__class__.__name__}')
        return self.__class__(*self.image(other.elems))
    compose = __mul__

    def __pow__(self, exp : int) -> 'Permutation':
        if exp == 0:
            return self.__class__(*self.elems) # identity composition, return copy of self to avoid mutation
        elif exp < 0:
            return self.inverse.__pow__(abs(exp))
        else:
            return reduce(mul, [self for _ in range(exp)]) # TODO : see if there's a more efficient way to do this
        
    @property
    def order(self) -> int:
        '''Smallest power of a permutation which generates the identity permutation'''
        return lcm(*(len(cycle) for cycle in self.to_cycles()))

    # MATRICES
    @staticmethod
    def _is_valid_permutation_matrix(matrix : np.ndarray) -> bool:
        return (                                  # A valid permutation matrix must be:
            (matrix.ndim == 2)                        # 1) 2-dimensional
            and (matrix.shape[0] == matrix.shape[1])  # 2) square
            and (np.in1d(matrix, [0 ,1])).all()       # 3) contain ONLY ones and zeroes
            and (matrix.sum(axis=0) == 1).all()       # 4) have exactly a single 1 in each column
            and (matrix.sum(axis=1) == 1).all()       # 4) have exactly a single 1 in each row
        )

    @classmethod
    def from_matrix(cls, matrix : np.ndarray) -> 'Permutation':
        '''Create a permutation from a permutation matrix'''
        if not cls._is_valid_permutation_matrix(matrix):
            raise ValueError('Provided matrix is not a valid permutation matrix')
        
        _, elems = np.nonzero(matrix.T)
        return cls(*elems)
        # _, elems = np.nonzero(matrix) # alternate, equally-valid implementation
        # return cls(*elems).inverse
    
    def to_matrix(self) -> np.ndarray:
        '''Obtain permutation matrix representation'''
        return np.eye(self.degree, dtype=int)[:, self.elems]
    
    @property
    def matrix(self) -> np.ndarray:
        return self.to_matrix()

    # WORD FORMS
    @classmethod
    def from_word(cls, word : str, delimiter : str='') -> 'Permutation':
        '''Create a permutation from a string with a total ordering'''
        return cls(*np.argsort(word.split(delimiter))).inverse # need to invert to produce expected result, as argsort gives the permutation which restores the identity

    def to_word(self) -> str:
        return ' '.join(str(i) for i in self.elems) # TODO : find best way to deal with multi-digit numbers

    # CYCLES
    def _cycle_decomposition(self) -> list[Cycle]:
        '''Return the disjoint cycles of a Permutation'''
        visited : list[bool] = [False] * self.degree
        cycles = []
        for elem in self.elems:
            cycle_elems = []
            while not visited[elem]:
                cycle_elems.append(elem)
                visited[elem] = True
                elem = self.elems[elem]
            if cycle_elems:
                cycles.append(Cycle(*cycle_elems))

        return cycles
    
    @classmethod
    def from_cycle(cls, cycle : Cycle, degree : Optional[int]=None) -> None: 
        '''Create a permutation from an ordering of elements in a single cycle. Missing elements are inferred to be fixed points'''
        max_elem = max(cycle) # determine value once to avoid multiple calls
        if degree is None:
            degree = max_elem
        elif (max_elem > degree - 1): # check if the degree is provided BUT is incompatible with the natural number set provided
            raise ValueError(f'Cannot construct a permutation of degree {degree} for a set containing at least {max_elem + 1} elements')
        return cls(*(cycle.mapping.get(i, i) for i in range(degree + 1))) # insert mapped element if provided, or else one with the same index (i.e. a fixed point)

    @classmethod
    def from_cycles(cls, cycles : list[Cycle]) -> 'Permutation':
        '''Stitch together Permutation from disjoint cycle representation'''
        degree = Cycle.max_elem_in_cycles(cycles) + 1
        return reduce(
            mul,
            (cls.from_cycle(cycle, degree=degree) for cycle in cycles[::-1]), # apply cycles right-to-left
            cls.identity(degree) # initialize with identity permutation
        ) 
    
    def as_cycle(self) -> Cycle:
        '''Reinterpret permutation a defining a cyclic order'''
        return Cycle(*self.elems)

    def to_cycles(self, canonicalize : bool=True) -> list[Cycle]:
        '''Produce cycle decompositon of permutation. NOT TO BE CONFUSED WITH Permutation.as_cycles()!
        By default, returns in canonical order (i.e. each cycle is presented with the least element first)'''
        if not canonicalize:
            return self._cycle_decomposition()
        
        # implicit else for canonicalization of cycles
        return sorted( # in canonical cycle notation, cycles are:
            (cycle.starting_from_index(np.argmax(cycle)) for cycle in self._cycle_decomposition()), # 1) listed with the largest element first          
            key=lambda x : x[0]                                                            # 2) sorted in ascending order by the first element
        ) 

    @property
    def cycles(self) -> list[Cycle]:
        '''Cycle decompositon of permutation in the order that elements appear (i.e. NOT in canonical order)'''
        return self.to_cycles(canonicalize=False)

    @property
    def cycles_canonical(self) -> list[Cycle]:
        '''Cycle decompositon of permutation in canonical order (largest element) first in each cycle and)'''
        '''Get cycle decompositon of permutation is canonical order'''
        return self.to_cycles(canonicalize=True)
    standard_notation = standard_form = canonical_cycles = cycles_canonical # lots of aliases to avoid getting bogged down in notation

    @property
    def cycle_type(self) -> dict[int, int]:
        '''Returns a dict of cycle lengths and the number of cycle in the permutation with that length'''
        return Cycle.cycle_type(self.to_cycles())

    # INVERSIONS AND LEHMER CODES
    @classmethod
    def from_lehmer(cls, lehmer_code : Sequence[int]) -> 'Permutation':
        '''Stitch together permutation from an inversion vector'''
        elems = [elem for elem in reversed(lehmer_code)] # create reversed copy to simplify iteration indexing avoid modifying original lehmer code
        for i, elem in enumerate(elems):
            for j, right_elem in enumerate(elems[:i]):
                if right_elem >= elem:
                    elems[j] += 1
        return cls(*reversed(elems)) # reverse final result to recover permutation order

    def to_lehmer_code(self) -> list[int]:
        '''Construct the left lehmer code for a permutation
        Each position in the code gives the number of inversions to the right of the permutation element at the corresponding position'''
        vector = [elem for elem in self.elems] # initialize witha  copy of the positions
        for i, elem in enumerate(vector):
            for j, right_elem in enumerate(vector[i:], start=i): # start indices at i to preserve correct absolute position in list
                if right_elem > elem:
                    vector[j] -= 1
        return vector
    
    @property
    def lehmer_code(self) -> list[int]:
        return self.to_lehmer_code()
    inversion_vector = lehmer = lehmer_code # aliases for convenience

    @property
    def num_inversions(self) -> int:
        '''Get total number of "out-of-order" elements in the permutation'''
        return sum(self.lehmer)

    @property
    def sign(self) -> int:
        '''The parity of the number of inversion of a permutation'''
        return -1 if (self.num_inversions % 2) else 1
    parity = sign

    @property
    def is_even(self) -> bool:
        return (self.sign == 1)

    @property
    def is_odd(self) -> bool:
        return (self.sign == -1)

    # ASCENTS AND DESCENTS
    @property
    def ascents(self) -> list[int]:
        '''The positions of elements which are followed by a greater element'''
        return [
            i for i, elem in enumerate(self.elems[:-1])
                if elem < self.elems[i + 1]
        ]
        
    @property
    def num_ascents(self) -> int:
        '''The number of elements which are followed by a greater element'''
        return len(self.ascents)

    @property
    def descents(self) -> list[int]:
        '''The positions of elements which are followed by a lesser element'''
        return [
            i for i, elem in enumerate(self.elems[:-1])
                if elem > self.elems[i + 1]
        ]
        
    @property
    def num_descents(self) -> int:
        '''The number of elements which are followed by a lesser element'''
        return len(self.descents)

    @property
    def support(self) -> list[int]:
        '''The positions of elemnts which do not map to themselves'''
        return [i for i, elem in enumerate(self.elems) if elem != i]
    
    @property
    def support_size(self) -> int:
        '''The number of elements in the permutation which do not map to themselves'''
        return len(self.support)
    
    # PERMUTATION GROUPS
    @staticmethod
    def cycle_index(perm_group : Iterable['Permutation'], variable : str='x', add_sym : str=' + ', mul_sym : str='*', exp_sym : str='^') -> str:
        '''Construct a string of the cycle index of a permutation group'''
        group_size : int = 0 # this allows passing of Iterable permutation groups and avoid needing to know group size a priori
        cycle_type_counts = defaultdict(int)
        for perm in perm_group:
            cycle_index_str = Cycle.cycle_index_sym(perm.cycle_type, variable, mul_sym, exp_sym)
            cycle_type_counts[cycle_index_str] += 1
            group_size += 1

        return f'1/{group_size}*({add_sym.join(f"{count}{mul_sym}{cycle_index_str}" for cycle_index_str, count in cycle_type_counts.items())})'

    @classmethod
    def symmetric_group(cls, order : int) -> Generator['Permutation', None, None]:
        '''Generate all permutations of a given order'''
        for elems in permutations(range(order), order):
            yield cls(*elems)

    @classmethod
    def alternating_group(cls, order : int) -> Generator['Permutation', None, None]:
        '''Generate all even permutations of a given order'''
        raise NotImplemented

    @classmethod
    def cyclic_group(cls, order : int) -> Generator['Permutation', None, None]:
        '''Generate all cyclic permutations of a given order'''
        raise NotImplemented

    @classmethod
    def dihedral_group(cls, order : int) -> Generator['Permutation', None, None]:
        '''Generate all permutations of vertex-numbered regular polygon'''
        raise NotImplemented