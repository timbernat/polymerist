'''Utilities for checking and enforcing module dependencies within code'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Callable, ParamSpec, TypeVar

Params = ParamSpec('Params')
ReturnType = TypeVar('ReturnType')
TCall = Callable[Params, ReturnType] # generic function of callable class

import importlib
from functools import wraps


def module_installed(module_name : str) -> bool:
    '''
    Check whether a module of the given name is present on the system
    
    Parameters
    ----------
    module_name : str
        The name of the module, as it would occur in an import statement
        Do not support direct passing of module objects to avoid circularity 
        (i.e. no reason to check if a module is present if one has already imported it elsewhere)
    
    Returns
    -------
    module_found : bool
        Whether or not the module was found to be installed in the current working environment
    '''
    try:
        package = importlib.import_module(module_name)
    except ModuleNotFoundError:
        return False
    else:
        return True
    
def modules_installed(*module_names : list[str]) -> bool:
    '''
    Check whether one or more modules are all present
    Will only return true if ALL specified modules are found
    
    Parameters
    ----------
    module_names : *str
        Any number of module names, passed as a comma-separated sequence of strings
        
    Returns
    -------
    all_modules_found : bool
        Whether or not all modules were found to be installed in the current working environment
    '''
    return all(module_installed(module_name) for module_name in module_names)

def requires_modules(*required_module_names : list[str]) -> Callable[[TCall[..., ReturnType]], TCall[..., ReturnType]]:
    '''
    Decorator which enforces optional module dependencies prior to function execution
    
    Parameters
    ----------
    module_names : *str
        Any number of module names, passed as a comma-separated sequence of strings
        
    Raises
    ------
    ImportError : Exception
        Raised if any of the specified packages is not found to be installed
        Exception message will indicate the name of the specific package found missing
    '''
    def decorator(func) -> TCall[..., ReturnType]:
        @wraps(func)
        def req_wrapper(*args : Params.args, **kwargs : Params.kwargs) -> ReturnType:
            for module_name in required_module_names:
                if not module_installed(module_name):
                    raise ImportError(f'No installation found for module "{module_name}"')
            else:
                return func(*args, **kwargs)
            
        return req_wrapper
    return decorator