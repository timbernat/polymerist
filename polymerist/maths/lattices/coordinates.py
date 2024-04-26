'''Representation of vectors of coordinates and elementary distance geometry operations'''

from typing import Optional, TypeVar, Union
from ...genutils.typetools.numpytypes import Shape, N, M, Num

import numpy as np
from itertools import product as cartesian_product


# CLASSES
class Coordinates:
    '''Represention class for encapsulating M-vectors of N-dimensional coordinates, with basic coordinate operations'''
    def __init__(self, points : Optional[np.ndarray[Shape[M, N], Num]]=None) -> None:
        if points is None:
            points = np.array([[]])

        assert(points.ndim == 2) # only accepts vector of coordinate
        self.points = points
        self.n_dims = points.shape[-1] # the number of dimensions of a point in the vector

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self.n_dims}-dimensional vector of {self.n_points} points)'

    # COORDINATE VECTOR DIMENSIONS AND EXTREMAL POINTS
    @property
    def n_points(self) -> int:
        '''The number of points currently stored in the coordinates vector'''
        return self.points.shape[0] # len(self.points)

    @property
    def dimensions(self) -> np.ndarray[Shape[N], int]:
        '''The side lengths of the bounding box along each dimension'''
        return self.points.ptp(axis=0)
    dims = sidelens = sidelengths = dimensions

    @property
    def minimum(self) -> np.ndarray[Shape[N], Num]:
        '''The bounding box vertex with the smallest coordinates in each dimension'''
        return self.points.min(axis=0)
    min = lower = smallest = minimum

    @property
    def maximum(self) -> np.ndarray[Shape[N], Num]:
        '''The bounding box vertex with the largest coordinates in each dimension'''
        return self.points.max(axis=0)
    max = upper = largest = maximum

    @property
    def extrema(self) -> np.ndarray[Shape[2, N], Num]:
        '''A 2xN array of the minimal and maximal bounding box vertices'''
        return np.vstack([self.minimum, self.maximum])

    # POINT OPERATIONS
    def __call__(self, index : int) -> np.ndarray[Shape[N], Num]:
        '''Retrieve the point at the given index'''
        return self.points[index]

    def _point_is_compat(self, point : np.ndarray[Shape[N], Num]) -> bool:
        '''Internal method to check whether an incoming point is compatible with the dimensions of the current array'''
        return (
            isinstance(point, np.ndarray)
            and (point.ndim == 1)
            and (point.size == self.n_dims)
            # TODO: check compatible types? (may be too inflexible, categorical.Numeric currently doesn't recognize numpy-specific types)
        )
    
    def validate_point(self, point : np.ndarray[Shape[N], Num]) -> None:
        '''Check if a point is compatible with the coordinates, or else raise Exception with reason why not'''
        if not self._point_is_valid(point):
            raise ValueError(f'Incompatible point {point}')

    @property
    def dists_to_point(self, point : np.ndarray[Shape[N], Num], norm_order : Optional[int]=None) -> np.ndarray[Shape[M], Num]:
        '''The distance between each point in a coordinate array and a single arbitrary point'''
        assert(self._point_is_compat(point))
        return np.linalg.norm(self.points - point, ord=norm_order, axis=1)

    # MEASURES OF CENTRALITY
    def weighted_centroid(self, weights : Optional[np.ndarray[Shape[M], Num]]=None) -> np.ndarray[Shape[N], Num]:
        '''The average (center-of-mass) coordinate of a vector of coordinates, with optional array of weights for each coordinates'''
        if weights is None:
            weights = np.ones(self.n_points, dtype=self.points.dtype)
        
        # prechecks to ensure the weights vector is of compatible size
        assert(weights.ndim == 1)
        assert(weights.size == self.points.shape[0])
        weights = weights.reshape((self.points.shape[0], 1))
        
        return (self.points * weights).mean(axis=0)

    @property
    def centroid(self) -> np.ndarray[Shape[N], Num]:
        '''The unweighted average coordinate of the vector'''
        return self.weighted_centroid() # weighed centroid with default unit weights
    center_of_mass = COM = centroid

    @property
    def dists_to_centroid(self, norm_order : Optional[int]=None, weights : Optional[np.ndarray[Shape[M], Num]]=None) -> np.ndarray[Shape[M], Num]:
        '''The distance of each coordinate in an array of coordinates to the coordinates' centroid'''
        return self.dists_to_point(point=self.weighted_centroid(weights=weights), norm_order=norm_order)
    effective_radii = eff_rad = dists_to_centroid

    # POINT ORDERINGS
    @property
    def lex_ordered_idxs(self) -> np.ndarray[Shape[M, N], int]: # TOSELF: "int" is the correct type annotation here (NOT Num)
        '''Returns a vector of the position that each point in self.points occupies when ordered lexicographically'''
        return np.lexsort(self.points.T)

    @property
    def lex_ordered_points(self) -> np.ndarray[Shape[M, N], Num]:
        '''Return copy of the points in the lattice in lexicographic order'''
        return self.points[self.lex_ordered_idxs]

    def lex_order_points(self) -> None:
        '''Sort points in the lattice in lexicographic order'''
        self.points = self.lex_ordered_points

    def randomize_points(self) -> None:
        '''Place the points in the lattice in a random order'''
        np.random.shuffle(self.points)

    # LATTICE TRANSFORMATIONS
    def translate(self, displacement : np.ndarray[Shape[N], Num]) -> None:
        '''Apply affine shift (translation only) to all points'''
        self.validate_point(displacement)
        self.points += displacement # TODO: use explicit broadcast here to reduce ambiguity?

    def linear_transformation(self, matrix : np.ndarray[Shape[N,N], float], as_coords : bool=False) -> Union[np.ndarray[Shape[M, N], float], 'Coordinates']:
        '''Accepts an NxN matrix (where N is the dimension of the lattice), returns a linearly-transformed copy of the coordinate points'''
        assert(matrix.shape == (self.n_dims, self.n_dims))
        transformed_points = self.points @ matrix.T # NOTE: need to right-multiply and transpose, since ROWS of self.points need to be tranformed 

        if as_coords:
            return self.__class__(transformed_points)
        return transformed_points

    def affine_transformation(self, matrix : np.ndarray[Shape[N,N], float], as_coords : bool=False) -> Union[np.ndarray[Shape[M, N], float], 'Coordinates']: # TOSELF: typehint on input matrix should be of shape N+1, N+1
        '''Accepts an NxN matrix (where N is the dimension of the lattice), returns an affinely-transformed copy of the coordinate points'''
        assert(matrix.shape == (self.n_dims + 1, self.n_dims + 1))
        aug_points = np.concatenate([self.points, np.ones((self.n_points, 1), dtype=int)], axis=1) # augment points vectors with extra columns of ones
        aug_transformed = aug_points @ matrix.T
        transformed_points = aug_transformed[: , :self.n_dims] / aug_transformed[:, self.n_dims, None] # downcast augmented transformed points from homogeneous coordinates, normalizing by projective part

        if as_coords:
            return self.__class__(transformed_points)
        return transformed_points


class BoundingBox(Coordinates):
    '''For representing a minimum bounding box around an M-coordinate vector of N-dimensional points'''
    def __init__(self, coords : Union[Coordinates, np.ndarray[Shape[M, N], Num]]) -> None:
        if isinstance(coords, np.ndarray):
            coords = Coordinates(coords) # allow passing of Coordinates-like classes
        points = np.array([vertex for vertex in cartesian_product(*self.extrema.T)])
        
        super().__init__(points)

    @property
    def vertices(self) -> np.ndarray[Shape[M, N], Num]: # TOSELF : first axis of returned array is actually of size 2**N (haven't implemented typing yet)
        '''Alias of the "points" attribute to distinguish use via syntax'''
        return self.points
    
    def surrounds(self, coords : Union[Coordinates, np.ndarray[Shape[M, N], Num]], strict : bool=False) -> np.ndarray[Shape[M, N], bool]:
        '''Boolean mask of whether the coordinates in a point vector lies within the bounding box'''
        if isinstance(coords, np.ndarray):
            coords = Coordinates(coords)
        assert(coords.n_dims == self.n_dims)
        less_funct = np.less if strict else np.less_equal # set "<" vs "<=" check by strict flag

        return less_funct(self.lower, coords.points) & less_funct(coords.points, self.upper)