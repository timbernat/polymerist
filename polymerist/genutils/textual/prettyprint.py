'''For generating human-readable string representations of other Python objects'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Any

from textwrap import indent
from enum import StrEnum


class Justification(StrEnum):
    '''For specifying string justification'''
    LEFT   = '<'
    CENTER = '^'
    RIGHT  = '>'
Just = Justification # alias for the lazy or hurried

def procrustean_string(
        string : str,
        length : int,
        padding : str=' ',
        just : Justification=Justification.LEFT,
    ) -> int:
    '''Takes a string and a target length and returns a new string which begins
    with the same characters as the original string but is clamped to the target length,
    truncating or padding if the original string is too long or short, respectively
    
    Parameters
    ----------
    string : str
        The string to stretch or cut
    length : int
        The target number of characters in the final string
    padding : str, default=" "
        A single character which shold be used as padding 
        when strings are too short, by default just a space
        MUST BE EXACTLY ONE CHARACTER!
    just : Justification, default=Justification.LEFT
        Enum specifier of how to justify a padded string
        Options are Justification.LEFT, Justification.CENTER, or Justification.RIGHT  
        
    Returns
    -------
    fmt_str : str
        A string which begins with the same characters as "string" but has
        precisely the specified length, with specified padding as specified
    '''
    if not (isinstance(length, int) and (length >= 0)):
        raise ValueError(f'Target string length must be a non-negative integer (not {length})')
    if not len(padding) == 1:
        raise IndexError(f'Padding string must contain exactly one character (passed "{padding}")')
        
    return f'{string[:length]:{padding}{just.value}{length}}'

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