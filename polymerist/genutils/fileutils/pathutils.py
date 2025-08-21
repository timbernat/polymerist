'''Utilities for editing, augmenting, and querying Paths'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Callable, Concatenate, ParamSpec, TypeVar, Union
T = TypeVar('T')
Params = ParamSpec('Params')

from inspect import signature
from functools import wraps
 
from pathlib import Path

from ..decorators.meta import extend_to_methods
from ..decorators.signatures import modify_param_annotation_by_index
    

# PATH CONVERSION FUNCTIONS (FOR CHANGING BETWEEN TYPES)
def aspath(path : Union[str, Path]) -> Path:
	'''Allow functions which expect Paths to also accept strings'''
	if not isinstance(path, Path):
		path = Path(path)
	return path

def asstrpath(strpath : Union[str, Path]) -> str:
	'''Allow functions which expect strings paths to also accept Paths'''
	if not isinstance(strpath, str):
		strpath = str(strpath)
	return strpath

@extend_to_methods
def allow_string_paths(funct : Callable[[Concatenate[Path, Params]], T]) -> Callable[[Concatenate[Union[Path, str], Params]], T]:
    '''Modifies a function which expects a Path as its first argument to also accept string-paths'''
    # TODO : add assertion that the wrapped function has at least one arg AND that the first arg is of the desired (limited) type
    old_sig = signature(funct) # lookup old type signature

    @wraps(funct) # for preserving docstring and type annotations / signatures
    def str_path_wrapper(flex_path : Union[str, Path], *args : Params.args, **kwargs : Params.kwargs) -> T:
        '''First converts stringy paths into normal Paths, then executes the original function'''
        return funct(aspath(flex_path), *args, **kwargs)

    # MODIFY SIGNATURE OF PATH-LIKE FIRST ARGUMENT TO MATCH NEW TYPE FLEXIBILITY
    str_path_wrapper.__signature__ = modify_param_annotation_by_index(
        old_sig,
        index=0, # modify signature of first argument to reflect new type flexibility
        new_type=Union[Path, str]
    )
    return str_path_wrapper

@extend_to_methods
def allow_pathlib_paths(funct : Callable[[Concatenate[str, Params]], T]) -> Callable[[Concatenate[Union[Path, str], Params]], T]:
    '''Modifies a function which expects a string path as its first argument to also accept canonical pathlib Paths'''
    # TODO : add assertion that the wrapped function has at least one arg AND that the first arg is of the desired (limited) type
    old_sig = signature(funct) # lookup old type signature

    @wraps(funct) # for preserving docstring and type annotations / signatures
    def str_path_wrapper(flex_path : Union[str, Path], *args : Params.args, **kwargs : Params.kwargs) -> T:
        '''First converts normal Paths into stringy paths, then executes the original function'''
        return funct(asstrpath(flex_path), *args, **kwargs)

    # MODIFY SIGNATURE OF PATH-LIKE FIRST ARGUMENT TO MATCH NEW TYPE FLEXIBILITY
    str_path_wrapper.__signature__ = modify_param_annotation_by_index(
        old_sig,
        index=0, # modify signature of first argument to reflect new type flexibility
        new_type=Union[Path, str]
    )
    return str_path_wrapper


# EMPTINESS CHECKS
@allow_string_paths
def is_empty_dir(dirpath : Path) -> bool:
    '''Check if a directory contains no files'''
    if not dirpath.is_dir():
         raise NotADirectoryError(f'dirpath must point to directory, not to file "{dirpath}"')
    return not any(dirpath.iterdir()) # can't use "len" for generators : TODO : make this more efficient (i.e. iteration-based) for large directories

@allow_string_paths
def is_empty_file(filepath : Path) -> bool:
    '''Check if a file contains no data'''
    if filepath.is_dir():
        raise IsADirectoryError(f'filepath must point to file, not to directory "{filepath}"')
    # NOTE: not checking file existence here, as calling stat() will already do this check (and raise appropriate error)

    return filepath.stat().st_size == 0

@allow_string_paths
def is_empty(path : Path) -> bool:
    '''Flexibly check whether a path is "empty"
    If path point to a file, returns whether the file contains data
    If path points to a directory, returns whether the directory contains any files (empty or otherwise)
    '''
    if path.is_dir():
        return is_empty_dir(path)
    elif path.is_file():
        return is_empty_file(path)
    else:
        raise FileNotFoundError(f'Path "{path}" does not exist')
    
# PATH CREATION FUNCTIONS
@allow_string_paths
def assemble_path(directory : Path, prefix : str, extension : str, postfix : str='') -> Path:
    '''Combine output, naming, descriptive, and filetype info to generate a complete Path'''
    if extension[0] == '.':
        extension = extension[1:] # remove leading dots if included
    path_name = f'{prefix}{"_" if postfix else ""}{postfix}.{extension}'

    return directory / path_name

# PATH PROPERTY FUNCTIONS (DON'T MODIFY ANYTHING)
def _dotless(extension : str) -> str:
    '''Separate the dot from a SINGLE extension file suffix. Returns the original suffix if not dot is present'''
    return extension.split('.')[-1] 

@allow_string_paths
def dotless(path : Path) -> str:
    '''Separate the dot from file path. Returns the original suffix if not dot is present'''
    return _dotless(path.suffix)

# PATH MODIFICATION FUNCTIONS (CHANGING THE STRUCTURE OF A PATH OBJECT)
@allow_string_paths
def default_suffix(path : Path, suffix : str) -> Path:
    '''Asserts that a path has a suffix, appending a specified default suffix if none exists'''
    if not path.suffix:
        path = path.with_name(f'{path.stem}.{suffix}') # ensure charge params path has correct extension

    return path

@allow_string_paths
def prepend_parent(path : Path, new_parent : Path) -> Path:
    '''Prepends a parent tree to an existing path'''
    return new_parent / path

@allow_string_paths
def detach_parent(path : Path, old_parent : Path) -> Path:
    '''Cuts off a parent tree from an existing path'''
    return path.relative_to(old_parent)

@allow_string_paths
def exchange_parent(path : Path, old_parent : Path, new_parent : Path) -> Path:
    '''Exchanges the parent tree of a path for another parent tree'''
    return prepend_parent(path=detach_parent(path, old_parent), new_parent=new_parent)

@allow_string_paths
def local_rename(path : Path, new_name : str) -> Path:
    '''Performs file rename relative to the parent directory (NOT the cwd)'''
    return path.rename(path.with_name(new_name))

@allow_string_paths
def local_restem(path : Path, new_stem : str) -> Path:
    '''Performs file rename relative to the parent directory (NOT the cwd), preserving the extension of the original file'''
    return path.rename(path.with_stem(new_stem))
