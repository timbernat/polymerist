'''Specialized cases of general positional numbering systems which are more common in usage'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Callable, Generator, Iterable
from .positional import PositionalNumberingSystem


def hypergeometric_ratios(funct : Callable[[int], int]) -> Generator[int, None, None]:
    '''Generates ratios between successive natural number terms according to a provided function'''
    i = 0
    while True:
        yield funct(i)
        i += 1


class FixedRadixNumberSystem(PositionalNumberingSystem):
    '''Positional numbering system with a single fixed radix'''
    def __init__(self, radix : int=10) -> None:
        self._radix = radix

    @property
    def radix(self) -> int:
        return self._radix
    base = radix

    @property
    def radices(self) -> Iterable[int]: # override radices getter with generator
        return hypergeometric_ratios(lambda i : self.radix)

class FactorialNumberSystem(PositionalNumberingSystem):
    '''For representing factoradic numbers, useful in enmerating permutations via Lehmer codes'''
    def __init__(self) -> None:
        pass

    @property
    def radices(self) -> Iterable[int]:
        return hypergeometric_ratios(lambda i : i + 1)
