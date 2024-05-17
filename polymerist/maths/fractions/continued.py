'''Representations and computation methods for continued fractions and ration approximations to real numbers'''

from typing import Union, Generator, Iterator, TypeAlias, Type
Shape : TypeAlias = tuple
Real : TypeAlias = Union[float, int]

import numpy as np


# CONSTANT PARAMETERS SHARED AMONGST MANY FUNCTIONS BELOW
INT_TYPE : Type = np.int64
EPS = 1E-8
TOL = 1E-6

# CONTINUED FRACTION CALCULATIONS
def real_to_continued_fraction_coeffs(x : Real, eps : float=EPS) -> Generator[int, None, None]:
    '''Euclidean algorithm on a real number to produce the integral continued fraction coefficients'''
    while True:
        n, rem = divmod(x, 1)
        if isinstance(n, float):
            assert(n.is_integer())
        # yield int(n)
        yield INT_TYPE(n)
        
        if abs(rem - 0) < eps:
            break # TOSELF : can't make part of while condition (without redundant pre-code) due to initial Euclidean division step
        x = 1 / rem

def continuant_matrix(a : int) -> np.ndarray[tuple[2, 2], int]:
    '''Homographic matrix for a evaluating the next continuant from a continued fraction coefficient'''
    return np.array([
        [a, 1],
        [1, 0]
    ], dtype=INT_TYPE)

def continued_fraction_to_continuants(coeffs : Iterator[int]) -> Generator[tuple[int, int], None, None]:
    '''Fold a sequence of continued fraction coefficients into successive continuants (raional approximations)'''
    extraction_vector = np.array([1,0], dtype=INT_TYPE).reshape(-1, 1)
    M = np.eye(2, dtype=INT_TYPE)
    for c in coeffs:
        M = M @ continuant_matrix(c)
        col = M @ extraction_vector
        yield tuple(col[:, 0])

# RATIONAL APPROXIMATIONS FROM CONTINUED FRACTIONS
def rational_approxes(x : Real, tol : float=TOL, eps : float=EPS) -> Generator[tuple[int, int], None, None]:
    '''Unfold a real number into its continued fraction representation, then generate successive continuants of it'''
    for p, q in continued_fraction_to_continuants(real_to_continued_fraction_coeffs(x, eps=eps)):
        yield p, q
        if abs(p/q - x) < tol:
            break

def best_rational_approx(x : Real, tol : float=TOL, eps : float=EPS) -> tuple[int, int]:
    '''Provide a rational approximation to a value with the smallest denominator that is within some tolerance'''
    *_, best = rational_approxes(x, tol=tol, eps=eps)
    return best
