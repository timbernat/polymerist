'''For graphical enumeration, i.e counting numbers of distinct graphs with various constraints'''

from functools import lru_cache
from .core import binomial_coeff


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
