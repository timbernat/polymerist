'''Tools for matrix decomposition'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import numpy as np
from ...genutils.typetools.numpytypes import Shape, DType, N, M


def diagonalize(matrix : np.ndarray) -> tuple[np.ndarray]:
    '''Diagonalize a matrix into it's eigenbasis. Return rotation and diagonal matrices P, D, and P^-1''' 
    eivals, eivecs = np.linalg.eig(matrix) # perform eigendecomposition to obtain principle components
    P = eivecs # eigenvector matrix is conveniently already in the form we desire
    D = np.diag(eivals)
    Pinv = np.linalg.inv(eivecs) # compute inverse matrix to rotate back

    return P, D, Pinv 

def inv_left(matrix : np.ndarray[Shape[M, N], DType]) -> np.ndarray[Shape[N, M], DType]:
    '''Return the left-inverse of an arbitrary (not necessarily square) matrix'''
    return np.linalg.inv(matrix.T @ matrix) @ matrix.T

def inv_right(matrix : np.ndarray[Shape[M, N], DType]) -> np.ndarray[Shape[N, M], DType]:
    '''Return the right-inverse of an arbitrary (not necessarily square) matrix'''
    return matrix.T @ np.linalg.inv(matrix @ matrix.T)