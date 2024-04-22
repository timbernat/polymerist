'''Type-hinting for Union-bound categories of types'''

from typing import Any, Callable, Container, Iterable, Sequence, Type, TypeAlias, Union
from numpy import number as np_number
from numpy.typing import ArrayLike

import builtins
from inspect import isclass
from .parametric import U


# TYPECHECKING FOR CUSTOM UNION TYPE ALIASES
def _union_member_factory(union : U) -> Callable[[Any], bool]:
    '''Factory for making Uion-membership-checking functions'''
    def isinunion(var : Any) -> bool:
        '''Check if an object is a member of a Union class'''
        return isinstance(var, union.__args__)
    return isinunion

Numeric : TypeAlias = Union[int, float, complex, np_number]
StringLike : TypeAlias = Union[str, bytes, bytearray]
JSONSerializable : TypeAlias = Union[str, bool, int, float, tuple, list, dict]

_TYPE_UNIONS : tuple[Type] = (Numeric, StringLike, JSONSerializable)
for klass in _TYPE_UNIONS:
    globals[f'is{klass.__name__.lower()}'] = _union_member_factory(klass) # register to module-level scope


# REGISTRIES OF ALL BUILTIN TYPES WITH GIVEN BASE CLASS BEHAVIOR
BUILTIN_TYPES : dict[str, Type] = {}
_BUILTIN_BASES_TO_CHECK : tuple[Type] = (Exception, Callable, Container, Iterable, Sequence)
for klass in _BUILTIN_BASES_TO_CHECK:
    class_reg : dict[str, Type] = {}
    for builtin_name in dir(builtins):
        builtin_obj = getattr(builtins, builtin_name)
        if isclass(builtin_obj):
            BUILTIN_TYPES[builtin_obj.__name__] = builtin_obj
            if issubclass(builtin_obj, klass):
                class_reg[builtin_obj.__name__] = builtin_obj
    globals()[f'BUILTIN_{klass.__name__.upper()}S']  = class_reg # register to module-level scope