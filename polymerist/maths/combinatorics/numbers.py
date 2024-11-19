'''Utilities for calculating fundamental combinatorial numbers'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Iterable
from operator import mul
from math import factorial # not worth re-implementing here, the C-implementation is plenty fast
from functools import lru_cache, reduce


# Sequential utility functions
def alt_sign(n : int) -> int:
    '''Returns the value of -1^n (useful for alternating sums)'''
    # return 2*(n % 2) - 1 # arithmetic turns out to be slightly slower than explicit check
    return -1 if (n % 2) else 1 


# Binomial and multinomial coefficients
@lru_cache
def binomial_coeff(n : complex, k : int) -> float: # appears to be quicker than scipy.special.binom
    '''Calculates a generalized binomial coefficient
    Counts the number of distinct ways to make and unordered selection of k elemnts from a set of n items'''
    if (k < 0):
        return 0.0
    
    if (k > n//2) and isinstance(n, int):
        return binomial_coeff(n, n - k)

    bincoeff = 1.0
    for i in range(k):
        bincoeff *= (n - i) / (k - i)

    return bincoeff
choose = binomial_coeff # alias for convenience

def multiset_coeff(n : int, k : int) -> float:
    '''Calculates a multiset coefficient
    Counts the number of ways to form a multiset of cardinality n from a support of size k'''
    return binomial_coeff(n + k -1, k)

def multinomial_coeff(multips : Iterable[int]) -> float:
    '''Calculate the multinomial coefficient of a set of multiplicities
    Counts the number of distinct permutations of a multiset with given multiplicites'''
    num = factorial(sum(multips))
    denom = reduce(mul, (factorial(i) for i in multips)) 

    return num / denom

def multinomial_coeff_native(multips : Iterable[int]) -> float:
    '''Calculate the multinomial coefficient of a set of multiplicities
    Counts the number of distinct permutations of a multiset with given multiplicites

    Native Python implementation of multinomial_coeff() which is slightly slower but requires no imports'''
    multicoeff, total = 1, 0 # initialize with multiplicate and additive identities, respectively
    for elem in multips:     # use expansion as product of binomials
        total += elem
        multicoeff *= binomial_coeff(total, elem)

    return multicoeff

def catalan(n : int) -> float:
    '''Calculates the n-th Catalan number
    Counts dozens of different things, including the number of monotonic 2D random walks
    and the number of valid (i.e. fully closed) groupings of n pairs of parentheses'''
    return binomial_coeff(2*n, n) / (n + 1)


# Counting partitions and cycles
def pentagonal(n : int) -> float:
    '''Calculates the n-th pentagonal number
    Counts the number of points in n "concentric" pentagon figures which share a corner vertex

    Related to the generating function for numbers of partitions'''
    # return binomial_coeff(n, 1) + 3*binomial_coeff(n, 2)
    return (3*n**2 - n) / 2

@lru_cache
def stirling_second(n : int, k : int) -> float:
    '''Calculates the Stirling number of the second kind S(n, k)
    Counts the number of ways to partition n objects into k non-empty pairwise-disjoint subsets

    Recurrence used as listed in NIST DLMF (https://dlmf.nist.gov/26.8)'''
    return sum(
        alt_sign(k - i) * binomial_coeff(k, i) * i**n
            for i in range(k + 1)
    ) / factorial(k)

@lru_cache
def stirling_first(n : int, k : int) -> float: # note : needs to be defined second, as the Stirling numbers of the second kinds are part of the definition here
    '''Calculates the Stirling number of the first kind c(n, k)
    Counts the number of ways to to decompose a permutation of n objects into k cycles

    Recurrence used as listed in Abramowitz ans Stegun (https://www.convertit.com/Go/ConvertIt/Reference/AMS55.ASP?Res=150&Page=825)'''
    if (n == 0) and (k == 0):
        return 1.0 # for some reason, the recurrence given does not give the correct boundary value

    return sum(
        alt_sign(j - k) * binomial_coeff(j - 1, k - 1) * binomial_coeff(2*n - k, j) * stirling_second(j - k, j - n)
            for j in range(n, 2*n - k + 1)
    )
    
@lru_cache
def bell(n : int) -> float:
    '''Calculates the n-th Bell number
    Counts the number of ways to partition n objects into any number of non-empty pairwise-disjoint subsets'''
    return sum(
        stirling_second(n, k)
            for k in range(n + 1)
    )

@lru_cache
def bernoulli(n : int) -> int: # NOTE: avoided implementing in terms of Stirling numbers of second kind due to catastophic cancellation of factorials
    '''Calculates the n-th Bernoulli number
    Related to the sums of powers of integers (via Faulhalber's formula)'''
    if n == 0:
        return 1
    return -sum( # NOTE: implemented recursively, as both explicit and Stirling number (2nd) implementations fail even for n < 20
        binomial_coeff(n + 1, k) * bernoulli(k)
            for k in range(n)
    ) / (n + 1)

# Graph enumeration
@lru_cache
def count_labelled_graphs(N : int) -> int:
    '''Gives number of all labelled graphs on N nodes'''
    return 2**binomial_coeff(N, 2)

@lru_cache
def count_connected_labelled_graphs(n : int) -> float:
    '''Gives array of number labelled graphs with one connected component on N nodes'''
    if n == 0:
        return 0.0 # TOSELF : should this be 1 instead? (Really more of a philosophical question, should the null graph be counted as having 1 component?)
    return count_labelled_graphs(n) - sum(
        binomial_coeff(n - 1, k - 1)*count_labelled_graphs(n - k)*count_connected_labelled_graphs(k)
            for k in range(1, n)
    )