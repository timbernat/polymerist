'''Additional type-hinting and typechecking not provided by builtins'''

from typing import Any, Callable, ParamSpec, TypeAlias, TypeVar, Union
import numpy.typing # NOTE : this import is necessary to use np.typing calls
import numpy as np


# GENERIC TYPE VARIABLES
T = TypeVar('T') # universal generic type

C = TypeVar('C') # generic type for a class
O = TypeVar('O') # generic type for an object passed to a function
F = TypeVar('F') # generic type for a function
U = TypeVar('U') # generic class from representing Literal Unions (since arg-free union is not supported)

P = ParamSpec('P') # for representing (preserved) input parameters
R = TypeVar('R')   # for representing generic return values

Args   = TypeVar('Args'  ) # generic type for arguments to a function
KWArgs = TypeVar('KWArgs') # generic type for keyword arguments to a function

# CATEGORICAL TYPEHINTS
def _union_member_factory(union : U) -> Callable[[Any], bool]:
    '''Factory for making Uion-membership-checking functions'''
    def isinunion(var : Any) -> bool:
        '''Check if an object is a member of a Union class'''
        return isinstance(var, union.__args__)
    return isinunion

Numeric = Union[int, float, complex, np.number]
ArrayLike = Union[list, tuple, np.ndarray, np.typing.ArrayLike] # exclude dicts and sets, as they introduce complicated behaviors
JSONSerializable = Union[str, bool, int, float, tuple, list, dict] 

isnumeric = _union_member_factory(Numeric)
isarraylike = _union_member_factory(ArrayLike)
isjsonserializable = _union_member_factory(JSONSerializable)

# NUMPY-SPECIFIC TYPEHINTS
N = TypeVar('N', bound=int) # for typing arbitrary array dimension
M = TypeVar('M', bound=int) # for typing arbitrary array dimension
DType = TypeVar('DType', bound=np.generic)
Shape : TypeAlias = tuple

Coords2D : TypeAlias = np.ndarray[Shape[M, 2], DType] # an Mx2 array of M 3D coordinates
Coords3D : TypeAlias = np.ndarray[Shape[M, 3], DType] # an Mx3 array of M 3D coordinates
CoordsND : TypeAlias = np.ndarray[Shape[M, N], DType] # an MxN array of M ND coordinates

Vector2D : TypeAlias = np.ndarray[Shape[2], DType] # a 1x2 array of M 3D coordinates
Vector3D : TypeAlias = np.ndarray[Shape[3], DType] # a 1x3 array of M 3D coordinates
VectorND : TypeAlias = np.ndarray[Shape[N], DType] # a 1xN array of M ND coordinates