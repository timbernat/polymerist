'''Custom linear algebra utilities not already supplied by a numeric processing library (like numpy)'''

import numpy as np


def diagonalize(matrix : np.ndarray) -> tuple[np.ndarray]:
    '''Diagonalize a matrix into it's eigenbasis. Return rotation and diagonal matrices P, D, and P^-1''' 
    eivals, eivecs = np.linalg.eig(matrix) # perform eigendecomposition to obtain principle components
    P = eivecs # eigenvector matrix is conveniently already in the form we desire
    D = np.diag(eivals)
    Pinv = np.linalg.inv(eivecs) # compute inverse matrix to rotate back

    return P, D, Pinv 