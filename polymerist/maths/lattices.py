'''Utilities for generating periodic unit lattices'''

from typing import Iterable

import numpy as np
from itertools import product as cartesian_product

def generate_int_lattice(*dims : Iterable[int]) -> np.ndarray:
    '''Generate all N-D coordinates of points on a integer lattice with the sizes of all D dimensions given'''
    return np.fromiter(
        iter=cartesian_product(*[
            range(d)
                for d in dims
        ]),
        dtype=np.dtype((int, len(dims)))
    )