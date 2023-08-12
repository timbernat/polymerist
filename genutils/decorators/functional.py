'''Decorators for modifying functions'''

from typing import Callable, Optional, Union

from copy import deepcopy
from pathlib import Path

from ..typetools import O, T, Args, KWArgs
from ..fileutils.pathutils import aspath, asstrpath


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

def allow_string_paths(funct : Callable[[Path, Args, KWArgs], T]) -> Callable[[Union[Path, str], Args, KWArgs], T]:
    '''Modifies a function which expects a Path as its first argument to also accept string-paths'''
    def str_path_wrapper(flex_path : Union[str, Path], *args : Args, **kwargs : KWArgs) -> T:
        '''First converts stringy paths into normal Paths, then executes the original function'''
        return funct(aspath(flex_path), *args, **kwargs)
    return str_path_wrapper

def allow_pathlib_paths(funct : Callable[[str, Args, KWArgs], T]) -> Callable[[Union[Path, str], Args, KWArgs], T]:
    '''Modifies a function which expects a string path as its first argument to also accept canonical pathlib Paths'''
    def str_path_wrapper(flex_path : Union[str, Path], *args : Args, **kwargs : KWArgs) -> T:
        '''First converts normal Paths into stringy paths, then executes the original function'''
        return funct(asstrpath(flex_path), *args, **kwargs)
    return str_path_wrapper