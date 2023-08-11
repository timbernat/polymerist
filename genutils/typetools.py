'''Additional type-hinting and typechecking not provided by builtins'''

from typing import Any, TypeVar, Union
import numpy.typing # NOTE : this import is necessary to use np.typing calls
import numpy as np


# GENERIC TYPE VARIABLES
T = TypeVar('T') # universal generic type

C = TypeVar('C') # generic type for a class
O = TypeVar('O') # generic type for an object passed to a function
F = TypeVar('F') # generic type for a function

Args   = TypeVar('Args'  ) # generic type for arguments to a function
KWArgs = TypeVar('KWArgs') # generic type for keyword arguments to a function

# CATEGORICAL TYPEHINTS
Numeric = Union[int, float, complex, np.number]
def isnumeric(var : Any) -> bool:
	'''Check if a variable is numerical'''
	return isinstance(var, Numeric.__args__)

ArrayLike = Union[list, tuple, np.ndarray, np.typing.ArrayLike] # exclude dicts and sets, as they introduce complicated behaviors
def isarraylike(var : Any) -> bool:
	'''Check if a variable behaves as an indexable array'''
	return isinstance(var, ArrayLike.__args__)

JSONSerializable = Union[str, bool, int, float, tuple, list, dict] 