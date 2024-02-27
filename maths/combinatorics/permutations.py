'''Utilites for representing pure permutations, cycles, and permutation groups'''

from typing import Generator, Iterable, Optional, Sequence, TypeVar


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