'''Core tools for manipulating integer lattices in D-dimensions'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Iterable
from numbers import Number
from ...genutils.typetools.numpytypes import Shape, D, N

import numpy as np
from itertools import product as cartesian_product
from .coordinates import Coordinates, BoundingBox


def generate_int_lattice(*dims : Iterable[int]) -> np.ndarray[Shape[N, D], int]:
    '''Generate all D-D coordinates of points on a integer lattice with the sizes of all D dimensions given'''
    return np.fromiter(
        iter=cartesian_product(*[
            range(d)
                for d in dims
        ]),
        dtype=np.dtype((int, len(dims)))
    )

def nearest_int_coord_along_normal(point : np.ndarray[Shape[D], Number], normal : np.ndarray[Shape[D], Number]) -> np.ndarray[Shape[D], int]:
    '''
    Takes an D-dimensional control point and an D-dimensional normal vector
    Returns the integer-valued point nearest to the control point which lies in
    the normal direction relative to the control point
    '''
    # Check that inputs have vector-like shapes and compatible dimensions
    assert(point.ndim == 1)
    assert(normal.ndim == 1)
    assert(point.size == normal.size)

    min_int_bound_point = np.floor(point)
    max_int_bound_point = np.ceil( point)
    if np.isclose(min_int_bound_point, max_int_bound_point).all(): # edge case: if already on an integer-valued point, the fllor and ceiling will be identical
        int_point = point 
    else:
        int_bbox = BoundingBox(np.vstack([min_int_bound_point, max_int_bound_point]))  # generate bounding box from extremal positions
        dots = np.inner((int_bbox.vertices - point), normal) # take dot product between the normal and the direction vectors from the control point to the integer bounding points
        i = dots.argmax() # position of integer point in most similar direction to normal
        furthest_point, similarity = int_bbox.vertices[i], dots[i]

        if similarity <= 0:
            raise ValueError(f'Could not locate valid integral point in normal direction (closest found was {furthest_point} with dot product {similarity})')
        else:
            int_point = furthest_point 

    return int_point.astype(int)

# TODO : implement enumeration of integral points within an D-simplex

class CubicIntegerLattice(Coordinates[int]):
    '''For representing an n-dimensional integer lattice, consisting of all n-tuples of integers with values constrained by side lengths in each dimension'''
    def __init__(self, counts_along_dims : np.ndarray[Shape[D], int]) -> None: # TODO: implement more flexible input support (i.e. star-unpacking, listlikes, etc.)
        assert(counts_along_dims.ndim == 1)
        super().__init__(generate_int_lattice(*counts_along_dims))
        self.counts_along_dims = counts_along_dims

    def counts_along_dims_as_str(self, multip_char : str='x') -> str:
        '''Stringify the lattice sidelengths'''
        return multip_char.join(str(i) for i in self.counts_along_dims)

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self.n_dims}-dimensional, {self.counts_along_dims_as_str()})'

    # LATTICE DIMENSIONS
    @property
    def capacity(self) -> int: # referred to as "N" in typehints
        '''The maximum number of points that the lattice could contains'''
        return np.prod(self.counts_along_dims)

    @property
    def lex_ordered_weights(self) -> np.ndarray[Shape[D], int]:
        '''Vector of the number of points corresponding
        Can be viewed as a linear transformation between indices and point coordinates when coords are placed in lexicographic order'''
        return np.concatenate(([1], np.cumprod(self.counts_along_dims)[:-1]))

    # SUBLATTICE DECOMPOSITION
    @property
    def odd_even_idxs(self) -> tuple[np.ndarray[Shape[N], int], np.ndarray[Shape[N], int]]: # TOSELF: each subarray actually has length N/2 (+1 if capacity is odd), not sure how to typehint that though
        '''Return two vectors of indices, corresponding to points in the "odd" and "even" non-neighboring sublattices, respectively'''
        parity_vector = np.mod(self.points.sum(axis=1), 2) # remainder of sum of coordinates of each point; corresponds to the condition that a single step along any dimension should invert parity
        is_odd = parity_vector.astype(bool) # typecast as booleans to permit indexing (and make intent a bit clearer)

        return np.flatnonzero(is_odd), np.flatnonzero(~is_odd) # need to faltten to avoid inconvenient tuple wrapping

    @property
    def odd_idxs(self) -> np.ndarray[Shape[N], int]:
        '''Indices of the point in the in "odd" sublattice'''
        return self.odd_even_idxs[0]

    @property
    def even_idxs(self) -> np.ndarray[Shape[N], int]:
        '''Indices of the point in the in "even" sublattice'''
        return self.odd_even_idxs[1]

    @property
    def odd_sublattice(self) -> np.ndarray[Shape[N, D], int]:
        '''Returns points within the odd sublattice of the lattice points'''
        return self.points[self.odd_idxs]
    odd_points = odd_sublattice # alias for convenience

    @property
    def even_sublattice(self) -> np.ndarray[Shape[N, D], int]:
        '''Returns points within the even sublattice of the lattice points'''
        return self.points[self.even_idxs]
    even_points = even_sublattice # alias for convenience