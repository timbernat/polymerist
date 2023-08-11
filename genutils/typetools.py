'''Additional type-hinting and typechecking not provided by builtins'''

from typing import Any, Union
import numpy.typing # NOTE : this import is necessary to use np.typing calls
import numpy as np


Numeric = Union[int, float, complex, np.number]
def isnumeric(var : Any) -> bool:
	'''Check if a variable is numerical'''
	return isinstance(var, Numeric.__args__)

ArrayLike = Union[list, tuple, np.ndarray, np.typing.ArrayLike] # exclude dicts and sets, as they introduce complicated behaviors
def isarraylike(var : Any) -> bool:
	'''Check if a variable behaves as an indexable array'''
	return isinstance(var, ArrayLike.__args__)

JSONSerializable = Union[str, bool, int, float, tuple, list, dict] 