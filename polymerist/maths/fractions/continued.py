'''Representations and computation methods for continued fractions and ration approximations to real numbers'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Union, Generator, Iterable, Type, TypeAlias, TypeVar
from ...genutils.typetools.numpytypes import Shape, NPInt

import numpy as np

Real : TypeAlias = Union[float, int]
I = TypeVar('I', bound=Union[int, NPInt])


# CONSTANT DEFAULT PARAMETERS SHARED AMONGST MANY FUNCTIONS BELOW
DEFAULT_INT_TYPE : Type = np.int64
DEFAULT_EPS = 1E-8
DEFAULT_TOL = 1E-6

# CONTINUED FRACTIONS AND PARTIAL CONVERGENTS ("CONTINUANTS")
def continuant_matrix(a : int, int_type : I=DEFAULT_INT_TYPE) -> np.ndarray[Shape[2, 2], int]:
    '''Homographic matrix for a evaluating the next continuant from a continued fraction coefficient'''
    return np.array([
        [a, 1],
        [1, 0]
    ], dtype=int_type)

def _continuant_matrices_from_coeffs(coeffs : Iterable[int], int_type : I=DEFAULT_INT_TYPE) -> Generator[np.ndarray[Shape[2, 2], I], None, None]:
    '''Generate a sequence of continuant matrices representing adjacent steps in unravelling Euclidean division'''
    M = np.eye(2, dtype=int_type) # initialize first two pairs as [1, 0] and [0, 1]
    for c in coeffs:
        M = M @ continuant_matrix(c)
        yield M

def continued_fraction_to_continuants(coeffs : Iterable[int], int_type : I=DEFAULT_INT_TYPE) -> Generator[tuple[int, int], None, None]:
    '''Fold a sequence of continued fraction coefficients into successive continuants (rational approximations)'''
    for M in _continuant_matrices_from_coeffs(coeffs, int_type=int_type):
        yield M[:, 0]

# EUCLIDEAN ALGORITHM VARIANTS
def real_to_continued_fraction_coeffs(x : Real, eps : float=DEFAULT_EPS, int_type : I=DEFAULT_INT_TYPE) -> Generator[int, None, None]:
    '''Euclidean algorithm on a real number to produce the integral continued fraction coefficients'''
    while True:
        n, rem = divmod(x, 1)
        if isinstance(n, float):
            assert(n.is_integer())
        # yield int(n)
        yield int_type(n)
        
        if abs(rem - 0) < eps:
            break # TOSELF : can't make part of while condition (without redundant pre-code) due to initial Euclidean division step
        x = 1 / rem

def extended_euclidean_algorithm(a : int, b : int, int_type : I=DEFAULT_INT_TYPE) -> tuple[int, int, int]:
    '''Compute the Extended Euclidean Algorithm between two integers "a" and "b"
    Returns the greatest common divisor of a and b, along with a pair of Bezout coefficients"x" and "y" which satisfy a*x + b*y = gcd(a, b)'''
    quotients : list[int] = []
    ai, bi = a, b # make copies to avoid mutation and to have access to original values after iteration
    while bi:
        quot, rem = divmod(ai, bi)
        ai, bi = bi, rem
        quotients.append(-quot) # NOTE: negative sign is intentional here

    _gcd = ai # NOTE: underline is to avoid name clash with builtin math.gcd
    *_, bezout_coeffs = _continuant_matrices_from_coeffs(quotients, int_type=int_type) # NOTE: do have access to original a, b values in first column, but signs are annoying to keep track of
    y, x = bezout_coeffs[:, -1]
    assert((a*x + b*y) == _gcd) # verify that Bezout's identity is actually satisfied

    return _gcd, x, y
bezout_coeffs = euclidean_extended = extended_euclidean_algorithm

# RATIONAL APPROXIMATIONS FROM CONTINUED FRACTIONS
def rational_approxes(x : Real, tol : float=DEFAULT_TOL, eps : float=DEFAULT_EPS) -> Generator[tuple[int, int], None, None]:
    '''Unfold a real number into its continued fraction representation, then generate successive continuants of it'''
    for p, q in continued_fraction_to_continuants(real_to_continued_fraction_coeffs(x, eps=eps)):
        yield p, q
        if abs(p/q - x) < tol:
            break

def best_rational_approx(x : Real, tol : float=DEFAULT_TOL, eps : float=DEFAULT_EPS) -> tuple[int, int]:
    '''Provide a rational approximation to a value with the smallest denominator that is within some tolerance'''
    *_, best = rational_approxes(x, tol=tol, eps=eps)
    return best
