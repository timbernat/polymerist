'''Utilities for examining prime numbers and integer factorizations'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import TypeAlias

from math import sqrt
from operator import mul

from functools import reduce
from collections import defaultdict


Factorization : TypeAlias = dict[int, int] # dictionary of factors and their exponents

def is_prime(n : int) -> bool: # faster (but slightly less pretty) implementation
    '''Check if an integer is prime'''
    if n < 2:
        return False

    i = 2
    while (i * i) <= n:
        if (n % i) == 0:
            return False
        i += 1
    return True # must be prime if no smaller divisors are found

def is_prime_alt(n : int) -> bool: # more readable and Pythonic but slightly slower (same asymptotic complexity however)
    '''Check if an integer is prime'''
    if n < 2:
        return False

    for i in range(2, int(sqrt(n)) + 1):
        if (n % i) == 0:
            return False
    else:
        return True # must be prime if no smaller divisors are found

def prime_factorization(n : int) -> Factorization:
    '''Computes prime factorization of an integer n. Returns factorization as a dict,
    where the keys are the prime factors and the values are their respective exponents'''
    factors = defaultdict(int)
    
    i = 2
    while (i * i) <= n:
        while (n % i) == 0:
            n /= i
            factors[i] += 1
        i += 1

    if n > 2:
        factors[int(n)] = 1 # add any leftover large prime factors after division
    return factors

def num_from_factorization(fac : Factorization) -> int:
    '''Reconstruct a number as a composition of its factors'''
    return reduce(mul, (base**exp for base, exp in fac.items()))