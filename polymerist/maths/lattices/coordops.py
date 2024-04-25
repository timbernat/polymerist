'''Distance geometry methods for arrays of coordinates'''

'''Utilities for generating periodic unit lattices'''

from typing import Optional, TypeVar
import numpy as np
from itertools import product as cartesian_product

from ...genutils.typetools.numpytypes import Shape, N, M, DType
from ...genutils.typetools.categorical import Numeric
T = TypeVar('T', bound=Numeric)


# COORDINATE-TO-POINT DISTANCES
def mean_coord(coords : np.ndarray[Shape[M, N], DType], weights : Optional[np.ndarray[Shape[M, 1], DType]]=None) -> np.ndarray[Shape[N], DType]:
    '''The average (center-of-mass) coordinate of a vector of coordinates, with optional array of weightsfor each coordinates'''
    if weights is None:
        weights = np.ones((coords.shape[0], 1), dtype=coords.dtype)
    assert(weights.size == coords.shape[0])
    weights = weights.reshape((coords.shape[0], 1))
    
    return (coords * weights).mean(axis=0)
center_of_mass = COM = mean_coord

def dists_to_point(coords : np.ndarray[Shape[M, N], DType], point : np.ndarray[Shape[N], DType], norm_order : Optional[int]=None) -> np.ndarray[Shape[M], DType]:
    '''The distance between each point in a coordinate array and a single arbitrary point'''
    return np.linalg.norm(coords - point, ord=norm_order, axis=1)

def dists_to_centroid(coords : np.ndarray[Shape[M, N], DType], norm_order : Optional[int]=None) -> np.ndarray[Shape[M], DType]:
    '''The distance of each coordinate in an array of coordinates to the coordinates' centroid'''
    return dists_to_point(coords, point=mean_coord(coords), norm_order=norm_order)

# BOUNDING BOXES
class BoundingBox:
    '''For representing a minimum bounding box around an M-coordinate vector of N-dimensional points'''
    def __init__(self, points : np.ndarray[Shape[M, N], T]) -> None:
        assert(points.ndim == 2) # only accepts vector of coordinate
        self.points = points
        self.n_points, self.n_dims = points.shape
        
    def clone(self) -> 'BoundingBox':
        '''Create a copy of the current bounding box'''
        return self.__class__(self.points)

    @property
    def dimensions(self) -> np.ndarray[Shape[N], int]:
        '''The side lengths of the bounding box along each dimension'''
        return self.points.ptp(axis=0)
    dims = sidelens = sidelengths = dimensions

    @property
    def minimum(self) -> np.ndarray[Shape[N], T]:
        '''The bounding box vertex with the smallest coordinates in each dimension'''
        return self.points.min(axis=0)
    min = lower = smallest = minimum

    @property
    def maximum(self) -> np.ndarray[Shape[N], T]:
        '''The bounding box vertex with the largest coordinates in each dimension'''
        return self.points.max(axis=0)
    max = upper = largest = maximum

    @property
    def extrema(self) -> np.ndarray[Shape[2, N], T]:
        '''A 2xN array of the minimal and maximal bounding box vertices'''
        return np.vstack([self.minimum, self.maximum])

    @property
    def vertices(self) -> np.ndarray[Shape[M, N], T]: # TOSELF : first axis of returned array is actually of size 2**N (haven't implemented typing yet)
        '''A full vector of the vertices of the bounding box'''
        return np.array([
            vertex
                for vertex in cartesian_product(*self.extrema.T)
        ])
    
    def surrounds(self, coords : np.ndarray[Shape[M, N], T], strict : bool=False) -> np.ndarray[Shape[M, N], bool]:
        '''Boolean mask of whether the coordinates in a point vector lies within the bounding box'''
        assert(coords.ndim == 2)
        assert(coords.shape[1] == self.n_dims)

        less_funct = np.less if strict else np.less_equal # set "<" vs "<=" check by strict flag
        return less_funct(self.lower, coords) & less_funct(coords, self.upper)

# NORMAL VECTORS
def nearest_int_coord_along_normal(point : np.ndarray[Shape[N], Numeric], normal : np.ndarray[Shape[N], Numeric]) -> np.ndarray[Shape[N], int]:
    '''
    Takes an N-dimensional control point and an N-dimensional normal vector
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
        dots = np.inner((int_bbox.extrema - point), normal) # take dot product between the normal and the direction vectors from the control point to the integer bounding points
        i = dots.argmax() # position of integer point in most similar direction to normal
        furthest_point, similarity = int_bbox.extrema[i], dots[i]

        if similarity <= 0:
            raise ValueError(f'Could not locate valid integral point in normal direction (closest found was {furthest_point} with dot product {similarity})')
        else:
            int_point = furthest_point 

    return int_point.astype(int)