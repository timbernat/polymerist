'''For checking whether object are valid Python modules and packages, and if so for gathering info from within them'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Optional, Union
from types import ModuleType
from pathlib import Path

from importlib.resources import (
    Package,
    files as get_package_path
)
from importlib.resources._common import get_package, from_package, resolve


# CHECKING PACKAGE AND MODULE STATUS
def is_module(module : Package) -> bool:
    '''Determine whether a given Package-like (i.e. str or ModuleType) is a valid Python module
    This will return True for packages, bottom-level modules (i.e. *.py) and Python scripts'''
    try:
        resolve(module)
        return True
    except ModuleNotFoundError:
        return False
    
def is_package(package : Package) -> bool:
    '''Determine whether a given Package-like (i.e. str or ModuleType) is a valid Python package'''
    try:
        get_package(package)
        return True
    except (ModuleNotFoundError, TypeError):
        return False

# EXTRACTING MODULE NAMING INFO
def flexible_module_pass(module : Union[str, Path, ModuleType]) -> ModuleType: # TODO: extend this to decorator
    '''Flexible interface for supplying a ModuleType object as an argument
    Allows for passing a name (either module name or string path), Path location, or a module proper'''
    if isinstance(module, (str, ModuleType)):
        return resolve(module)
    elif isinstance(module, Path):
        raise NotImplementedError
    else:
        raise TypeError(f'Cannot interpret object of type "{type(module).__name__}" as a module')

# TODO : find way to get depth of submodule in toplevel ("number of dots" before standalone name)
def module_parts(module : Union[str, ModuleType]) -> tuple[Optional[str], str]:
    '''Takes a module (as its name or as ModuleType) and returns its parent package name and relative module name'''
    module = resolve(module)
    module_name = module.__spec__.name
    parent_package_name, _, module_stem = module_name.rpartition('.') # split on rightmost dot separator
    if not parent_package_name:
        parent_package_name = None

    return parent_package_name, module_stem

def module_stem(module : Union[str, ModuleType]) -> tuple[Optional[str], str]:
    '''Takes a module (as its name or as ModuleType) and returns its relative module name'''
    return module_parts(module)[-1]

def relative_module_name(module : ModuleType, relative_to : Optional[ModuleType]=None, remove_leading_dot : bool=True) -> str:
    '''Gets the name of a module relative to another (presumably toplevel) module
    If the given module is not in the path of the toplevel module, will simply return as module.__name__'''
    rel_mod_name = module.__spec__.name
    if relative_to is not None:
        toplevel_prefix = relative_to.__spec__.name
        if remove_leading_dot:
            toplevel_prefix += '.' # append dot to prefix to remove it later
        rel_mod_name = rel_mod_name.removeprefix(toplevel_prefix)

    return rel_mod_name

# FETCHING RESOURCES FROM PATHS WITHIN PACKAGES
def get_resource_path_within_package(relative_path : Union[str, Path], package : Package) -> Path:
    '''Get the Path to a resource (i.e. either a directory or a file) which lives within a Python package'''
    package_path : Path = get_package_path(package) # will also implicitly check that the provided package exists as a module
    resource_path = package_path / relative_path    # concat to Path here means string inputs for relative_path are valid without explicit conversion

    if not resource_path.exists(): # if this block is reached, it means "package" is a real module and resource path is DEFINED relative to package's path, so the below message is valid
        raise ValueError(f'{resolve(package).__name__} contains no resource "{relative_path}"')
    
    return resource_path

def get_dir_path_within_package(relative_path : Union[str, Path], package : Package) -> Path:
    '''Get the Path to a directory which lives within a Python package'''
    dir_path : Path = get_resource_path_within_package(package=package, relative_path=relative_path) # performs all check associated with getting the resource
    
    if not dir_path.is_dir():
        raise NotADirectoryError(f'{resolve(package).__name__} contains "{dir_path}", but it is not a directory')
    
    return dir_path

def get_file_path_within_package(relative_path : Union[str, Path], package : Package) -> Path:
    '''Get the Path to a (non-directory) file which lives within a Python package'''
    file_path : Path = get_resource_path_within_package(package=package, relative_path=relative_path) # performs all check associated with getting the resource
    
    if not file_path.is_file():
        raise FileNotFoundError(f'{resolve(package).__name__} contains no file "{file_path}"')
    
    return file_path