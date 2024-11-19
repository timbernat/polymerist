'''For generating human-readable string representations of other Python objects'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Any
from textwrap import indent


def dict_to_indented_str(dict_to_stringify : dict[Any, Any], level_delimiter : str='\t', line_sep : str='\n') -> str:
    '''Generate a pretty-printable string from a (possibly nested) dictionary,
    with each level of nesting indicated by "level_delimiter"'''
    text = []
    for key, value in dict_to_stringify.items():
        if isinstance(value, dict):
            text.append(f'{key!r}')
            text.append(indent(dict_to_indented_str(value), level_delimiter)) # recursive call for nested dicts
        else:
            text.append(f'{key!r} : {value!r}') # call repr methods

    return line_sep.join(text)