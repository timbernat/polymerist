'''Type-hinting for Union-bound categories of types'''

from typing import Any, Callable, Literal, Type, Union
from numpy import number as np_number
from numpy.typing import ArrayLike

from .parametric import U


def _union_member_factory(union : U) -> Callable[[Any], bool]:
    '''Factory for making Uion-membership-checking functions'''
    def isinunion(var : Any) -> bool:
        '''Check if an object is a member of a Union class'''
        return isinstance(var, union.__args__)
    return isinunion

Numeric = Union[int, float, complex, np_number]
JSONSerializable = Union[str, bool, int, float, tuple, list, dict] 

isnumeric = _union_member_factory(Numeric)
isarraylike = _union_member_factory(ArrayLike)
isjsonserializable = _union_member_factory(JSONSerializable)