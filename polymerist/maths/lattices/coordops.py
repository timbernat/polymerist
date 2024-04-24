'''Distance geometry methods for arrays of coordinates'''

'''Utilities for generating periodic unit lattices'''

from typing import Optional
import numpy as np
from itertools import product as cartesian_product

from ...genutils.typetools.numpytypes import Shape, N, M, DType


def mean_coord(coords : np.ndarray[Shape[M, N], DType]) -> np.ndarray[Shape[N], DType]:
    '''The average (center-of-mass) coordinate of a vector of coordinates'''
    return coords.mean(axis=0)
center_of_mass = COM = mean_coord

def dists_to_point(coords : np.ndarray[Shape[M, N], DType], point : np.ndarray[Shape[N], DType], norm_order : Optional[int]=None) -> np.ndarray[Shape[M], DType]:
    '''The distance between each point in a coordinate array and a single arbitrary point'''
    return np.linalg.norm(coords - point, ord=norm_order, axis=1)

def dists_to_centroid(coords : np.ndarray[Shape[M, N], DType], norm_order : Optional[int]=None) -> np.ndarray[Shape[M], DType]:
    '''The distance of each coordinate in an array of coordinates to the coordinates' centroid'''
    return dists_to_point(coords, point=mean_coord(coords), norm_order=norm_order)


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
