'''Distance geometry methods for arrays of coordinates'''

'''Utilities for generating periodic unit lattices'''

from typing import Optional
import numpy as np
from itertools import product as cartesian_product

from ...genutils.typetools.numpytypes import Shape, N, M, DType
from ...genutils.typetools.categorical import Numeric


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
def bounding_box_dims(coords : np.ndarray[Shape[M, N], DType]) -> np.ndarray[Shape[N], DType]:
    '''The N-dimensional tight bounding box side lengths around an array of coordinates'''
    return coords.ptp(axis=0)

def bounding_box_extrema(coords : np.ndarray[Shape[M, N], DType]) -> tuple[np.ndarray[Shape[N], DType], np.ndarray[Shape[N], DType]]:
    '''The minimal and maximal coordinates of the N-dimensional tight bounding box around an array of coordinate'''
    return coords.min(axis=0), coords.max(axis=0)

def bounding_box_points(coords : np.ndarray[Shape[M, N], DType]) -> np.ndarray[Shape[M, N], DType]: # TOSELF : first axis of returned array is actually of size 2**N (haven't implemented typing yet)
    '''The corner points the N-dimensional tight bounding box around an array of coordinate'''
    return np.array([
        corner_point
            for corner_point in cartesian_product(*np.vstack(bounding_box_extrema(coords)).T)
    ])

def coords_inside_bbox(coords : np.ndarray[Shape[M, N], DType], lower : np.ndarray[Shape[1, N], DType], upper : np.ndarray[Shape[1, N], DType], strict : bool=False) -> np.ndarray[Shape[M, N], bool]:
    '''Boolean mask of whether coordinates are within the boundaries of some bounding box
    With strict=True, points on the boundary are not considered inside; with strict=False, points on the boundary are considered insde as well'''
    less_funct = np.less if strict else np.less_equal # set "<" vs "<=" check by strict flag
    return less_funct(lower, coords) & less_funct(coords, upper)

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
        int_bound_extrema = np.vstack([min_int_bound_point, max_int_bound_point]) # package extremal pointsinto single array
        int_bound_points = bounding_box_points(int_bound_extrema) # expand extremal points into complete N-dimensional integer box around the control point
        dots = np.inner((int_bound_points - point), normal) # take dot product between the normal and the direction vectors from the control point to the integer bounding points
        i = dots.argmax() # position of integer point in most similar direction to normal

        if dots[i] > 0:
            int_point = int_bound_points[i]
        else:
            raise ValueError(f'Could not locate valid integral point in normal direction (closest found was {int_bound_points[i]} with dot product {dots[i]})')
        pass

    return int_point.astype(int)