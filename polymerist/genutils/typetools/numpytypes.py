'''Type aliases specific to numpy arrays'''

from typing import TypeAlias, TypeVar
import numpy as np
import numpy.typing as npt
from .categorical import Numeric


D = TypeVar('D', bound=int) # for typing arbitrary array dimensionality
N = TypeVar('N', bound=int) # for typing arbitrary array dimension size
M = TypeVar('M', bound=int) # for typing arbitrary array dimension size
Num = TypeVar('Num', bound=Numeric) # for typing numerically-valued generics
DType = TypeVar('DType', bound=np.generic)
Shape : TypeAlias = tuple

Coords2D : TypeAlias = np.ndarray[Shape[M, 2], DType] # an Mx2 array of M 3D coordinates
Coords3D : TypeAlias = np.ndarray[Shape[M, 3], DType] # an Mx3 array of M 3D coordinates
CoordsND : TypeAlias = np.ndarray[Shape[M, N], DType] # an MxN array of M ND coordinates

Vector2D : TypeAlias = np.ndarray[Shape[2], DType] # a 1x2 array of M 3D coordinates
Vector3D : TypeAlias = np.ndarray[Shape[3], DType] # a 1x3 array of M 3D coordinates
VectorND : TypeAlias = np.ndarray[Shape[N], DType] # a 1xN array of M ND coordinates