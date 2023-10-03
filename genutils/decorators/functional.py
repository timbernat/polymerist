'''Decorators for modifying functions'''

from typing import Callable, Optional, Union

import inspect
from functools import wraps
from inspect import Parameter

from copy import deepcopy
from pathlib import Path

from .meta import extend_to_methods
from . import signatures
from ..typetools import O, T, Args, KWArgs
from ..fileutils.pathutils import aspath, asstrpath


@extend_to_methods
def optional_in_place(funct : Callable[[O, Args, KWArgs], None]) -> Callable[[O, Args, bool, KWArgs], Optional[O]]:
    '''Decorator function for allowing in-place (writeable) functions which modify object attributes
    to be not performed in-place (i.e. read-only), specified by a boolean flag'''
    # TODO : add assertion that the wrapped function has at least one arg AND that the first arg is of the desired (limited) type
    old_sig = inspect.signature(funct)
    
    @wraps(funct) # for preserving docstring and type annotations / signatures
    def in_place_wrapper(obj : O, *args : Args, in_place : bool=False, **kwargs : KWArgs) -> Optional[O]: # read-only by default
        '''If not in-place, create a clone on which the method is executed''' # NOTE : old_sig.bind screws up arg passing
        if in_place:
            funct(obj, *args, **kwargs) # default call to writeable method - implicitly returns None
        else:
            copy_obj = deepcopy(obj) # clone object to avoid modifying original
            funct(copy_obj, *args, **kwargs) 
            return copy_obj # return the new object
    
    # ADD IN-PLACE PARAMETER TO FUNCTION SIGNATURE
    new_sig = signatures.insert_parameter_at_index(
        old_sig,
        new_param=Parameter(
            name='in_place',
            default=False,
            annotation=bool,
            kind=Parameter.KEYWORD_ONLY
        ),
        index=signatures.get_index_after_positionals(old_sig)
    )

    # ANNOTATE MODIFED OBJECT RETURN TYPE AS Optional[<type>]
    mod_type = tuple(old_sig.parameters.values())[0].annotation # get annotation of the object being modified
    in_place_wrapper.__signature__ = new_sig.replace(return_annotation=Optional[mod_type]) # replace signature with appropriate modifications

    return in_place_wrapper

@extend_to_methods
def allow_string_paths(funct : Callable[[Path, Args, KWArgs], T]) -> Callable[[Union[Path, str], Args, KWArgs], T]:
    '''Modifies a function which expects a Path as its first argument to also accept string-paths'''
    # TODO : add assertion that the wrapped function has at least one arg AND that the first arg is of the desired (limited) type
    old_sig = inspect.signature(funct) # lookup old type signature

    @wraps(funct) # for preserving docstring and type annotations / signatures
    def str_path_wrapper(flex_path : Union[str, Path], *args : Args, **kwargs : KWArgs) -> T:
        '''First converts stringy paths into normal Paths, then executes the original function'''
        return funct(aspath(flex_path), *args, **kwargs)

    # MODIFY SIGNATURE OF PATH-LIKE FIRST ARGUMENT TO MATCH NEW TYPE FLEXIBILITY
    str_path_wrapper.__signature__ = signatures.modify_param_annotation_by_index(
        old_sig,
        index=0, # modify signature of first argument to reflect new type flexibility
        new_type=Union[Path, str]
    )

    return str_path_wrapper

@extend_to_methods
def allow_pathlib_paths(funct : Callable[[str, Args, KWArgs], T]) -> Callable[[Union[Path, str], Args, KWArgs], T]:
    '''Modifies a function which expects a string path as its first argument to also accept canonical pathlib Paths'''
    # TODO : add assertion that the wrapped function has at least one arg AND that the first arg is of the desired (limited) type
    old_sig = inspect.signature(funct) # lookup old type signature

    @wraps(funct) # for preserving docstring and type annotations / signatures
    def str_path_wrapper(flex_path : Union[str, Path], *args : Args, **kwargs : KWArgs) -> T:
        '''First converts normal Paths into stringy paths, then executes the original function'''
        return funct(asstrpath(flex_path), *args, **kwargs)

    # MODIFY SIGNATURE OF PATH-LIKE FIRST ARGUMENT TO MATCH NEW TYPE FLEXIBILITY
    str_path_wrapper.__signature__ = signatures.modify_param_annotation_by_index(
        old_sig,
        index=0, # modify signature of first argument to reflect new type flexibility
        new_type=Union[Path, str]
    )

    return str_path_wrapper