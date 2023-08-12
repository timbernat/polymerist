'''Additional type-hinting and typechecking not provided by builtins'''

from typing import Any, Callable, TypeVar, Union
import numpy.typing # NOTE : this import is necessary to use np.typing calls
import numpy as np


# GENERIC TYPE VARIABLES
T = TypeVar('T') # universal generic type

C = TypeVar('C') # generic type for a class
O = TypeVar('O') # generic type for an object passed to a function
F = TypeVar('F') # generic type for a function
U = TypeVar('U') # generic class from representing Literal Unions (since arg-free union is not supported)

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