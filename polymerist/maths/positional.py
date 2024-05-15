'''Conversion tools for representing positive integers in fixed and mixed radix positional bases'''

from typing import Any, Callable, Generator, Iterable, Sequence, Union
from math import inf


class PositionalNumberingSystem:
    '''For representing positive integers in standard and mixed-radix positional numbering systems'''
    def __init__(self, radices : Sequence[int]) -> None:
        '''Radices should be passed in ascending order of significance (i.e. least-significant bases first)'''
        self._radices = radices # add check for integers
     
    # descriptor for place values access (read-only)
    @property
    def radices(self) -> Sequence[int]:
        return self._radices
    bases = radices
    # NOTE: specifically opted not to include an iter() conversion property, as this would reset a generator on each call inside the loop

    @property
    def places(self) -> Generator[int, None, None]:
        '''Generates the values of each position in a general mixed-radix positional notation'''
        p = 1 # TOSELF: should also yield 1?
        for radix in self.radices:
            yield p
            p *= radix

    # implementations of general base conversions
    def int_to_digits_iter(self, n : int) -> Generator[int, None, None]:
        '''Convert a non-negative integer to its digit representation under the specified radices
        Yields digits in ascending order of significance (i.e. starting from least significant digit)'''
        if n < 0:
            raise ValueError('Cannot generate digits of a negative integer')

        radices = iter(self.radices) # convert to iterable to generator-like parsing of sequences
        while n:
            try:
                radix = next(radices)
            except StopIteration: # if radices have run out, yield residual as "infinitieth" digit
                yield n
                break

            n, digit = divmod(n, radix)
            if digit < 0: # account for downwards rounding when encountering negative radices
                n += 1
                digit -= radix

            yield digit

    def int_to_digits(self, n : int, as_str : bool=False) -> Union[str, list[int]]:
        '''Returns digits of number in descending order of significance (i.e. most significant digit first)'''
        digits = [digit for digit in self.int_to_digits_iter(n)]
        digits.reverse()
        
        if as_str:
            return ''.join(str(digit) for digit in digits)
        return digits

    def digits_to_int(self, digits : Sequence[int]) -> int:
        '''Convert a sequence of digits to an integer in the specified numbering system
        Digits should be passed in descending order (i.e. most-significant digit first)'''
        n : int = 0

        radices = iter(self.radices)
        for i, digit, place in enumerate(zip(reversed(digits), self.places), start=0):
            digit = int(digit) # allows for string passing

            try:
                radix = next(radices)
            except StopIteration:# radix is deliberately one position ahead in sequence, to allow "peek" ahead due to 
                radix = inf
            
            if not (0 <= digit <= radix):
                raise ValueError(f'Digit {digit} in position {i} exceeds max allowed digit value of {radix} for that place')
            n += digit*place
        
        return n 

    # dunders for concise calling
    def __call__(self, *args, **kwargs) -> list[int]:
        '''Syntactic sugar for digit-to-int conversion (expects a Sequence)'''
        if (len(args) == 1) and isinstance(args[0], Sequence):
            args = args[0]
        return self.digits_to_int(args, **kwargs)
            
    def __getitem__(self, value : Any) -> int:
        '''Syntactic sugar for int-to-digit conversion'''
        if isinstance(value, int):
            return self.int_to_digits(value, as_str=True)
        else:
            raise TypeError


# USEFUL NUMBER SYSTEMS WITH PARTICULAR SPECIALIZED RADICES
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
    

# initialization of some common bases
Factoradic = FactorialNumberSystem()

commmon_bases : dict[str, int] = {
    'binary'      : 2,
    'ternary'     : 3,
    'seximal'     : 6,
    'octal'       : 8,
    'decimal'     : 10,
    'duodecimal'  : 12,
    'hexadecimal' : 16,
}
for base_name, base in commmon_bases.items():
    base_sys = FixedRadixNumberSystem(base)
    globals()[base_name.capitalize()]     = base_sys # systematic name
    globals()[f'base{base}'.capitalize()] = base_sys # common name

    globals()[f'nega{base_name}'.capitalize()] = FixedRadixNumberSystem(-base) # also include negative base for kicks