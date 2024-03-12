'''Tools for manipulating and pretty-printing text from files and string-like objects'''

from typing import Any, Callable, Optional

import re
from pathlib import Path
from textwrap import indent

from .fileutils.extensions import FileTypeError


# CASE CONVERSION
def snake_case_to_camel_case(varname : str) -> str:
    '''Convert a name from Snake Case to Camel Case
    E.g. name_of_a_thing -> NameOfAThing'''
    return ''.join(word.capitalize() for word in varname.split('_'))

def camel_case_to_snake_case(varname : str) -> str:
    '''Convert a name from Camel Case to Snake Case
    E.g. NameOfAThing -> name_of_a_thing'''
    cap_idxs = [i for i, char in enumerate(varname) if char.isupper()]
    return '_'.join(
        varname[i_start:i_end].lower()
            for i_start, i_end in zip(cap_idxs, cap_idxs[1:]+[None])
    ) 


# STRING INTERPOLATION
def insert_into_text_periodic(text : str, period : int, insertion : str='\n') -> str:
    '''Takes a string of text and another "insertion" string and inserts it throughout the text every <period> characters'''
    return insertion.join(text[i:i+period] for i in range(0, len(text), period))

def insert_into_text_periodic_re(text : str, period : int, insertion : str='\n') -> str:
    '''Takes a string of text and another "insertion" string and inserts it throughout the text every <period> characters
    Same as insert_into_text_periodic(), but implemented with regular expressions (allows for more complicated logical extensions)'''
    SPACE_RE = re.compile(f'(?s)(.{{{period}}})') # double curly braces escape the f-string syntax (to use as literals in regex quanitifer)
    return re.sub(SPACE_RE, f'\\1{insertion}', text)


# ORDINALS
def ordinal_suffix_from_int(n : int) -> str:
    '''Produce the appropriate word suffix for an integer in sequential order
    E.g 1 -> "st" as in "first", 17 -> "th" as in "seventeenth, etc.'''
    n = int(abs(n)) # will tolerate negative values and floats which look like ints
    dsuff = {
        1 : 'st',
        2 : 'nd',
        3 : 'rd'
    }
    if (1 <= (i := n % 10) <= 3) and not (11 <= (n % 100) <= 13): # check for special case of 1, 2, and 3 EXCEPT when occuring in 11, 12, or, 13
        return dsuff[i]
    return 'th' # everything else ends in "th"

def ordinal_suffix_from_int_alt(n : int) -> str:
    '''Produce the appropriate word suffix for an integer in sequential order
    E.g 1 -> "st" as in "first", 17 -> "th" as in "seventeenth, etc.
    
    An alternative, slightly-slower but aesthetically-pleasing implementation''' # taken from SO https://stackoverflow.com/questions/9647202/ordinal-numbers-replacement/ 
    n = int(abs(n)) # will tolerate negative values and floats which look like ints
    if (n % 100) in (11, 12, 13): # annoying special cases for first few tens values:
        return 'th'
    return ['th', 'st', 'nd', 'rd', 'th'][min(n % 10, 4)] # min and extra "th" at end give sneaky way of defaulting to 'th'

def ordinal_from_int(n : int) -> str:
    '''Produce the word representation of an integer sequential order
    E.g 1 -> "1st", 17 -> "seventeenth", 33 -> "33rd", etc.'''
    return str(n) + ordinal_suffix_from_int(n)


# FORMATTING OTHER TYPES INTO STRINGS
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


# TEXT SEARCHING/EDITING
def filter_text_by_condition(in_text_path : Path, condition : Callable[[str], bool], out_text_path : Optional[Path]=None, postfix : str='filtered', inclusive : bool=True, return_filtered_path : bool=False) -> Optional[Path]:
    '''Create a copy of a text-based file containing only the lines which match to a given boolean condition
    
    If no explicit output path is given, will create an output file in the same directory as the source file
    with the same name plus "postfix" tacked on. Can optionally return the path to the filtered file (else None)

    "Inclusive" kw governs whether to write lines which DO or DON'T meet the condition'''
    if out_text_path is None:
        out_text_path = in_text_path.with_stem(f'{in_text_path.stem}{"_" if postfix else ""}{postfix}')

    if (out_text_path == in_text_path):
        raise PermissionError(f'Attempting to overwrite {in_text_path} with regex filter') # prevent write clash
    
    if (out_text_path.suffix != in_text_path.suffix):  # prevent file type conversion during transfer
        raise FileTypeError(f'Input and output file must have same extension (not {in_text_path.suffix} and {out_text_path.suffix})')

    with out_text_path.open('w') as outfile: 
        with in_text_path.open('r') as infile: # readfile is innermost in case error occurs during file read (caught by handler one level up)
            for line in infile:
                if (condition(line) == inclusive): # only write lines if (matching AND inclusive) OR (not matching AND exclusive)
                    outfile.write(line)

    if return_filtered_path:
        return out_text_path