'''Core tools for manipulating integer lattices in N-dimensions'''

from typing import Iterable
import numpy as np
from itertools import product as cartesian_product

from ...genutils.typetools.numpytypes import Shape, N, M


def generate_int_lattice(*dims : Iterable[int]) -> np.ndarray[Shape[M, N], int]:
    '''Generate all N-D coordinates of points on a integer lattice with the sizes of all D dimensions given'''
    return np.fromiter(
        iter=cartesian_product(*[
            range(d)
                for d in dims
        ]),
        dtype=np.dtype((int, len(dims)))
    )

class CubicIntegerLattice:
    '''For representing an n-dimensional integer lattice, consisting of all n-tuples of integers with values constrained by side lengths in each dimension'''
    def __init__(self, sidelens : np.ndarray[Shape[N], int]) -> None: # TODO: implement more flexible input support (i.e. star-unpacking, listlikes, etc.)
        assert(sidelens.ndim == 1)
        self.sidelens = sidelens # ordered vector of the number of points along each dimension
        self.points : np.ndarray[Shape[M, N], int] = generate_int_lattice(*self.sidelens)

    def sidelens_as_str(self, multip_char : str='x') -> str:
        '''Stringify the lattice sidelengths'''
        return multip_char.join(str(i) for i in self.sidelens)

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self.n_dims}-dimensional, {self.sidelens_as_str()})'

    # LATTICE DIMENSIONS
    @property
    def n_dims(self) -> int: # referred to as "N" in typehints
        '''The number of dimensions of the lattice'''
        return self.sidelens.size
    
    @property
    def capacity(self) -> int: # referred to as "M" in typehints
        '''The maximum number of points that the lattice could contains'''
        return np.prod(self.sidelens)

    @property
    def n_points(self) -> int: 
        '''The actual number of points currently contained in the lattice'''
        return self.points.shape[0]
    
    # LEXICOGRAPHIC ORDERING
    def __call__(self, index : int) -> np.ndarray[Shape[N], int]:
        '''Retrieve the point at the given index'''
        return self.points[index]

    @property
    def lex_ordered_weights(self) -> np.ndarray[Shape[N], int]:
        '''Vector of the number of points corresponding
        Can be viewed as a linear transformation between indices and point coordinates when coords are placed in lexicographic order'''
        return np.concatenate(([1], np.cumprod(self.sidelens)[:-1]))
    
    @property
    def lex_ordered_idxs(self) -> np.ndarray[Shape[M, N], int]:
        '''Returns a vector of the position that each point in self.points occupies when ordered lexicographically'''
        return np.lexsort(self.points.T)

    @property
    def lex_ordered_points(self) -> np.ndarray[Shape[M, N], int]:
        '''Return copy of the points in the lattice in lexicographic order'''
        return self.points[self.lex_ordered_idxs]

    # IN-PLACE POINT REORDERING
    def lex_order_points(self) -> None:
        '''Sort points in the lattice in lexicographic order'''
        self.points = self.lex_ordered_points

    def randomize_points(self) -> None:
        '''Place the points in the lattice in a random order'''
        np.random.shuffle(self.points)

    # SUBLATTICE DECOMPOSITION
    @property
    def odd_even_idxs(self) -> tuple[np.ndarray[Shape[M], int], np.ndarray[Shape[M], int]]: # TOSELF: each subarray actually has length M/2 (+1 if capacity is odd), not sure how to typehint that though
        '''Return two vectors of indices, corresponding to points in the "odd" and "even" non-neighboring sublattices, respectively'''
        parity_vector = np.mod(self.points.sum(axis=1), 2) # remainder of sum of coordinates of each point; corresponds to the condition that a single step along any dimension should invert parity
        is_odd = parity_vector.astype(bool) # typecast as booleans to permit indexing (and make intent a bit clearer)

        return np.where(is_odd), np.where(~is_odd)

    @property
    def odd_idxs(self) -> np.ndarray[Shape[M], int]:
        '''Indices of the point in the in "odd" sublattice'''
        return self.odd_even_idxs[0]

    @property
    def even_idxs(self) -> np.ndarray[Shape[M], int]:
        '''Indices of the point in the in "even" sublattice'''
        return self.odd_even_idxs[1]

    @property
    def odd_sublattice(self) -> np.ndarray[Shape[M, N], int]:
        '''Returns points within the odd sublattice of the lattice points'''
        return self.points[self.odd_idxs]
    odd_points = odd_sublattice # alias for convenience

    @property
    def even_sublattice(self) -> np.ndarray[Shape[M, N], int]:
        '''Returns points within the even sublattice of the lattice points'''
        return self.points[self.even_idxs]
    even_points = even_sublattice # alias for convenience

    # LATTICE TRANSFORMATIONS
    def linear_transformation(self, matrix : np.ndarray[Shape[N,N], float], periodic : bool=False) -> np.ndarray[Shape[M, N], float]:
        '''Accepts an NxN matrix (where N is the dimension of the lattice)
        Returns a linearly-transformed copy of the points currently in the lattice'''
        assert(matrix.shape == (self.n_dims, self.n_dims))
        return self.points @ matrix.T # NOTE: need to right-multiply and transpose, since ROWS of self.points need to be tranformed 

    def affine_transformation(self, matrix : np.ndarray[Shape[N,N], float], periodic : bool=False) -> np.ndarray[Shape[M, N], float]: # TOSELF: typehint on input matrix should be of shape N+1, N+1
        '''Accepts an (N+1)x(N+1) matrix (where N is the dimension of the lattice)
        Returns an affine-transformed copy of the points currently in the lattice'''
        assert(matrix.shape == (self.n_dims + 1, self.n_dims + 1))
        aug_points = np.concatenate([self.points, np.ones((self.n_points, 1), dtype=int)], axis=1) # augment points vectors with extra columns of ones
        aug_transformed = aug_points @ matrix.T

        return aug_transformed[: , :self.n_dims] / aug_transformed[:, self.n_dims, None] # downcast augmented transformed points from homogeneous coordinates, normalizing by projective part