'''Utilities for checking and enforcing module dependencies within code'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Callable, Optional, ParamSpec, TypeVar, Union

Params = ParamSpec('Params')
ReturnType = TypeVar('ReturnType')
TCall = Callable[Params, ReturnType] # generic function of callable class

# from importlib import import_module
from importlib.util import find_spec
from functools import wraps


class MissingPrerequisitePackage(Exception):
    '''Raised when a package dependency cannot be found and the user should be alerted with install instructions'''
    def __init__(self,
            importing_package_name : str,
            use_case : str,
            install_link : str,
            dependency_name : str,
            dependency_name_formal : Optional[str]=None
        ):
        if dependency_name_formal is None:
            dependency_name_formal = dependency_name
        
        message = f'''
        {use_case} require(s) {dependency_name_formal}, which was not found in the current environment
        Please install `{dependency_name}` by following the installation instructions at {install_link}
        Then try importing from "{importing_package_name}" again'''
        
        super().__init__(message)
        
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
    # try:
    #     package = import_module(module_name)
    # except ModuleNotFoundError:
    #     return False
    # else:
    #     return True
    
    try: # NOTE: opted for this implementation, as it never actually imports the package in question (faster and fewer side-effects)
        return find_spec(module_name) is not None
    except (ValueError, AttributeError, ModuleNotFoundError): # these could all be raised by a missing module
        return False
    
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

def requires_modules(
        *required_module_names : list[str],
        missing_module_error : Union[Exception, type[Exception]]=ImportError,
    ) -> Callable[[TCall[..., ReturnType]], TCall[..., ReturnType]]:
    '''
    Decorator which enforces optional module dependencies prior to function execution
    
    Parameters
    ----------
    module_names : *str
        Any number of module names, passed as a comma-separated sequence of strings
    missing_module_error : type[Exception], default ImportError
        The type of Exception to raise if a module is not found installed
        Defaults to ImportError
        
    Raises
    ------
    ImportError : Exception
        Raised if any of the specified packages is not found to be installed
        Exception message will indicate the name of the specific package found missing
    '''
    # meta-check to ensure type of raised Exception is valid
    if not isinstance(missing_module_error, Exception):
        if not (isinstance(missing_module_error, type) and issubclass(missing_module_error, Exception)):
            # DEV: this is potentially brittle, depending on how the specific Exception subtype is implemented?
            raise TypeError('Must pass either Exception instance or subtype to "missing_module_error') 
        
    def tailored_exception(module_name : str) -> Exception:
        '''Accessory function to generate targetted Exceptions based on the provided
        mssing_module_error value and the name of a module with no found installation'''
        if isinstance(missing_module_error, Exception):
            return missing_module_error
        
        if isinstance(missing_module_error, type):
           return missing_module_error(f'No installation found for module "{module_name}"')
    
    def decorator(func) -> TCall[..., ReturnType]:
        @wraps(func)
        def req_wrapper(*args : Params.args, **kwargs : Params.kwargs) -> ReturnType:
            for module_name in required_module_names:
                if not module_installed(module_name):
                    raise tailored_exception(module_name)
            else:
                return func(*args, **kwargs)
            
        return req_wrapper
    return decorator