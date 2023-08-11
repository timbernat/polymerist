'''Decorators for modifying functions'''

from typing import Callable, Optional, TypeVar
from copy import deepcopy


O = TypeVar('O') # generic type for an object passed to a function
F = TypeVar('F') # generic type for a function
Args   = TypeVar('Args'  ) # generic type for arguments to a function
KWArgs = TypeVar('KWArgs') # generic type for keyword arguments to a function

def optional_in_place(funct : Callable[[O, Args, KWArgs], None]) -> Callable[[O, Args, bool, KWArgs], Optional[O]]:
    '''Decorator function for allowing in-place (writeable) functions which modify object attributes
    to be not performed in-place (i.e. read-only), specified by a boolean flag'''
    def in_place_wrapper(obj : O, *args : Args, in_place : bool=False, **kwargs : KWArgs) -> Optional[O]: # read-only by default
        '''If not in-place, create a clone on which the method is executed'''
        if in_place:
            funct(obj, *args, **kwargs) # default call to writeable method - implicitly returns None
        else:
            copy_obj = deepcopy(obj) # clone object to avoid modifying original
            funct(copy_obj, *args, **kwargs) 
            return copy_obj # return the new object
    return in_place_wrapper