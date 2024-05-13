'''Conversion tools for representing the digits of numbers in various positional bases'''

from typing import Generator, Sequence
from math import log


# "standard" integer bases
def digitize(n : int, base : int=10, ascending : bool=False) -> Generator[int, None, None]:
    '''Yields the digits of an integer in a particular base (default base 10),
    in either ascending or non-ascending (i.e. descending, default) order of significance'''
    n = abs(n) # avoids infinite loop for negatives
    if n == 0: # avoid null return for zero values
        yield 0
    elif ascending:
        while n:
            n, d = divmod(n, base)
            yield d
    else:
        M = int(log(n, base)) # largest power of the base smaller than the actual number
        for i in range(M, -1, -1): # iterate from largest power of base 
            d, n = divmod(n, base**i)
            yield d


# factorial bases
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