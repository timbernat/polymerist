'''For dynamically inspecting and modifying attributes of Python objects'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Any, Optional, Union
import re


def compile_argfree_getable_attrs(obj : Any, getter_re : Union[str, re.Pattern]='.*', repl_str : Optional[str]=None) -> dict[str, Any]:
    '''Compile the values of all methods of an object which require no arguments other than perhaps the object itself (this EXCLUDES properties)
    Returns a dict whose keys are the names of the methods called and whose values are the return values of those object methods

    Can optionally filter the names of returned method using a regular expression, passed to "getter_re"
    Can also optionally replace the chosen regex with an arbitrary string (including the empty string), passed to "repl_str"
    
    Parameters
    ----------
    obj : Any
        Any object instance
    getter_re : str or re.Pattern (optional), default ".*"
        Optional regular expression to use for filtering down returned methods
        Only methods whose names match the target regex are returns
    repl_str : str (optional)
        If provided, will replace the 
        for example, repl_str="" can be used to delete the regex from returned method names

    Returns
    -------
    getable_dict : dict[str, Any]
        dict whose keys are the selected method names and whose values are the corresponding method returns
    '''
    getable_dict = {}
    for attr_name in dir(obj):
        if re.search(getter_re, attr_name): 
            try:
                attr_key = attr_name if (repl_str is None) else re.sub(getter_re, repl_str, attr_name)
                getable_dict[attr_key] = getattr(obj, attr_name)()
            except (TypeError, Exception): # TODO : find way to selectively intercept the Boost C++ wrapper ArgumentError
                pass
    return getable_dict