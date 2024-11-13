'''For dynamically inspecting and modifying attributes of Python objects'''

from typing import Any, Optional
import re


def compile_argfree_getable_attrs(obj : Any, getter_re : str='get', repl_str : Optional[str]=None) -> dict[str, Any]:
    '''Takes an object and returns a dict of the return values of all argument-free methods of the objects
    Looks for methods of the object whose names contain with "getter_re", and can replace this with the value of "repl_str" in the final dict output if provided'''
    getable_dict = {}
    for attr_name in dir(obj):
        if re.search(getter_re, attr_name): 
            try:
                attr_key = attr_name if (repl_str is None) else re.sub(getter_re, repl_str, attr_name)
                getable_dict[attr_key] = getattr(obj, attr_name)()
            except (TypeError, Exception): # TODO : find way to selectively intercept the Boost C++ wrapper ArgumentError
                pass
    return getable_dict