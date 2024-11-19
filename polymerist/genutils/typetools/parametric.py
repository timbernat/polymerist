'''Type aliases for Python callable input parameters'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import TypeVar, ParamSpec, _UnionGenericAlias

T = TypeVar('T') # universal generic type

C = TypeVar('C') # generic type for a class
O = TypeVar('O') # generic type for an object passed to a function
F = TypeVar('F') # generic type for a function
U = TypeVar('U', bound=_UnionGenericAlias) # generic class from representing Literal Unions (since arg-free union is not supported)

P = ParamSpec('P') # for representing (preserved) input parameters
R = TypeVar('R')   # for representing generic return values

# TOSELF: replace with typing.ParamSpecArgs/.ParamSpecKWArgs?
Args   = TypeVar('Args'  ) # generic type for arguments to a function
KWArgs = TypeVar('KWArgs') # generic type for keyword arguments to a function
