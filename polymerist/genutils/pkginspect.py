'''For checking whether object are valid Python modules and packages, and if so for gathering info from within them'''

from typing import Union
from pathlib import Path

from importlib.resources import (
    Package,
    files as get_package_path
)
from importlib.resources._common import get_package, from_package, resolve


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