'''Type-hinting for Union-bound categories of types'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Any, Callable, Container, Iterable, Sequence, Type, TypeAlias, Union
from numpy import (
    ndarray,
    number as np_number
)
from numpy.typing import ArrayLike

import builtins
from inspect import isclass
from .parametric import U


# TYPECHECKING FOR CUSTOM UNION TYPE ALIASES
def _union_member_factory(union : U, regname : str='Union') -> Callable[[Any], bool]:
    '''Factory for making Union-membership-checking functions'''
    def isinunion(var : Any) -> bool:
        return isinstance(var, union.__args__)
    
    isunion_function_name = f'is{regname.lower()}'
    isinunion.__name__ = isunion_function_name
    isinunion.__qualname__ = isunion_function_name
    isinunion.__doc__ = f'''Check if an object is an instance of the {regname} class'''

    return isinunion

_UNION_TYPES : dict[str, U] = { # registry of aliases type Unions (will be dynamically registered with isinstance-like checkers at module level)
    'StringLike'        : Union[str, bytes, bytearray],                         # classes which 
    'ListLike'          : Union[list, tuple, set, frozenset, bytes, bytearray], # classes which are Container-like and can be initialized from a list of values
    'JSONSerializable'  : Union[str, bool, int, float, tuple, list, dict],      # classes which can be written to a JSON file by the default JSONEncoder/Decoder
}
for union_name, union_type in _UNION_TYPES.items():
    globals()[union_name] = union_type
    _union_checker_func = _union_member_factory(union_type, union_name) # register to module-level scope
    globals()[_union_checker_func.__name__] = _union_checker_func # register to module-level scope


# REGISTRIES OF ALL BUILTIN TYPES WITH GIVEN BASE CLASS BEHAVIOR
BUILTIN_TYPES : dict[str, Type] = {}
_BUILTIN_BASES_TO_CHECK : tuple[Type] = (Exception, Callable, Container, Iterable, Sequence)
for _klass in _BUILTIN_BASES_TO_CHECK:
    class_reg : dict[str, Type] = {}
    for builtin_name in dir(builtins):
        builtin_obj = getattr(builtins, builtin_name)
        if isclass(builtin_obj):
            BUILTIN_TYPES[builtin_obj.__name__] = builtin_obj
            if issubclass(builtin_obj, _klass):
                class_reg[builtin_obj.__name__] = builtin_obj
    globals()[f'BUILTIN_{_klass.__name__.upper()}S']  = class_reg # register to module-level scope