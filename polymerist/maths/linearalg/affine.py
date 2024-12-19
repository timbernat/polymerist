'''Utilities to streamline creation of 4x4 affine transformation matrices of 3D linear transformations in homogeneous coordinates'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import numpy as np
import numpy.typing as npt
from typing import Annotated, Literal


# USEFUL CUSTOM TYPEHINTS
Array4x4 = Annotated[npt.NDArray[np.float64], Literal[4, 4]]
AffineMatrix = Array4x4 # alias for convenience

# RIGID TRANSFORMATIONS
def xyzTrans(x : float=0.0, y : float=0.0, z : float=0.0) -> AffineMatrix:
    '''Returns rigid affine matrix which performs an isometric translation by "x", "y", and "z" units
    Returns the Identity matrix in absence of passed arguments'''
    return np.array([
        [1, 0, 0, x],
        [0, 1, 0, y],
        [0, 0, 1, z],
        [0, 0, 0, 1],
    ], dtype='float64')


def xRot(angle_rad : float=0.0) -> AffineMatrix:
    '''Returns rigid affine matrix which rotates about the positive x-axis by "angle_rad" radians 
    Returns the Identity matrix in absence of passed arguments'''
    s = np.sin(angle_rad)
    c = np.cos(angle_rad)
    
    return np.array([
        [1, 0,  0, 0],
        [0, c, -s, 0],
        [0, s,  c, 0],
        [0, 0,  0, 1],
    ], dtype='float64')

def yRot(angle_rad : float=0.0) -> AffineMatrix:
    '''Returns rigid affine matrix which rotates about the positive y-axis by "angle_rad" radians 
    Returns the Identity matrix in absence of passed arguments'''
    s = np.sin(angle_rad)
    c = np.cos(angle_rad)
    
    return np.array([
        [c, 0, -s, 0],
        [0, 1,  0, 0],
        [s, 0,  c, 0],
        [0, 0,  0, 1],
    ], dtype='float64')

def zRot(angle_rad : float=0.0) -> AffineMatrix:
    '''Returns rigid affine matrix which rotates about the positive z-axis by "angle_rad" radians 
    Returns the Identity matrix in absence of passed arguments'''
    s = np.sin(angle_rad)
    c = np.cos(angle_rad)
    
    return np.array([
        [c, -s, 0, 0],
        [s,  c, 0, 0],
        [0,  0, 1, 0],
        [0,  0, 0, 1],
    ], dtype='float64')

def randRot(about_x : bool=True, about_y : bool=True, about_z : bool=True) -> AffineMatrix:
    '''Returns an affine matrix for a rotation by some random angle(s) about the x, y, and z axes (or any subset of those axes)'''
    rot_dir = { # concise way to encapsulate what rotations can be performed and whther to perform them
        xRot : about_x,
        yRot : about_y,
        zRot : about_z,
    }

    matrix = np.eye(4, dtype=float)
    for (rot_fn, should_rotate) in rot_dir.items():
        if should_rotate:
            matrix = rot_fn(2 * np.pi * np.random.rand()) @ matrix # generate random angle and multiply rotation into overall transform matrix (in order)

    return matrix


# SCALING AND SHEAR
def xyzScale(sx : float=1.0, sy : float=1.0, sz : float=1.0) -> AffineMatrix:
    '''Returns rigid affine matrix which scales basis vectors by factors of "sx", "sy", and "sz" units
    Returns the Identity matrix in absence of passed arguments'''
    return np.array([
        [sx,  0,  0, 0],
        [ 0, sy,  0, 0],
        [ 0,  0, sz, 0],
        [ 0,  0,  0, 1],
    ], dtype='float64')

# TODO : add shear matrix functions