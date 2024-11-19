'''For representing rational numbers, and more general ratios'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from dataclasses import dataclass
from typing import Any, Callable, ClassVar, TypeVar
N = TypeVar('N')

from numbers import Number
from math import gcd


# HELPER FUNCTIONS
def sgnmag(num : N) -> tuple[bool, N]:
    '''Returns the sign and magnitude of a numeric-like value'''
    return num < 0, abs(num)


# RATIO CLASSES
@dataclass(repr=False)
class Ratio:
    '''For representing fractional ratios between two objects'''
    num   : Any
    denom : Any

    # REPRESENTATION
    def __repr__(self) -> str:
        return f'{self.num}/{self.denom}'
    
    def to_latex(self) -> str:
        '''Return latex-compatible string which represent fraction'''
        return rf'\frac{{{self.num}}}{{{self.denom}}}'

    # RELATIONS
    @property
    def reciprocal(self) -> 'Ratio':
        '''Return the reciprocal of a ration'''
        return self.__class__(self.denom, self.num)

@dataclass(repr=False)
class Rational(Ratio):
    '''For representing ratios of integers'''
    num   : int
    denom : int

    # REDUCTION
    autoreduce : ClassVar[bool]=False
    
    def __post_init__(self) -> None:
        if self.__class__.autoreduce:
            self.reduce()

    def reduce(self) -> None:
        '''Reduce numerator and denominator by greatest common factor'''
        _gcd = gcd(self.num, self.denom)
        self.num=int(self.num / _gcd)
        self.denom=int(self.denom / _gcd)
    simplify = reduce # alias for convenience

    @property
    def reduced(self) -> 'Rational':
        '''Return reduced Rational equivalent to the current rational (does not modify in-place)'''
        new_rat = self.__class__(self.num, self.denom)
        new_rat.reduce()

        return new_rat
    simplifed = reduced # alias for convenience
    
    def as_proper(self) -> tuple[int, 'Rational']:
        '''Returns the integer and proper fractional component of a ratio'''
        integ, remain = divmod(self.num, self.denom)
        return integ, self.__class__(remain, self.denom)
    
    # ARITHMETIC
    def __add__(self, other : 'Rational') -> 'Rational':
        '''Sum of two Rationals'''
        return self.__class__(
            num=(self.num * other.denom) + (self.denom * other.num),
            denom=(self.denom * other.denom)
        )
    
    def __sub__(self, other : 'Rational') -> 'Rational':
        '''Difference of two Rationals'''
        return self.__class__(
            num=(self.num * other.denom) - (self.denom * other.num),
            denom=(self.denom * other.denom)
        )

    def __mul__(self, other : 'Rational') -> 'Rational':
        '''Product of two Rationals'''
        return self.__class__(
            num=self.num * other.num,
            denom=self.denom * other.denom
        )

    def __div__(self, other : 'Rational') -> 'Rational':
        '''Quotient of two Rationals'''
        return self.__class__(
            num=self.num * other.denom,
            denom=self.denom * other.num
        )
    
    def __pow__(self, power : float) -> 'Rational':
        '''Exponentiates a ratio'''
        return self.__class__(
            num=self.num**power,
            denom=self.denom**power
        )