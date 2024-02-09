'''Tools for manipulating and pretty-printing text from files and string-like objects'''

from typing import Any, Callable, Optional

from pathlib import Path
from textwrap import indent

from .fileutils.extensions import FileTypeError


# FORMATTING INTO STRINGS
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

# FILE TEXT EDITING
def filter_txt_by_condition(in_txt_path : Path, condition : Callable[[str], bool], out_txt_path : Optional[Path]=None, postfix : str='filtered', inclusive : bool=True, return_filtered_path : bool=False) -> Optional[Path]:
    '''Create a copy of a text-based file containing only the lines which match to a given boolean condition
    
    If no explicit output path is given, will create an output file in the same directory as the source file
    with the same name plus "postfix" tacked on. Can optionally return the path to the filtered file (else None)

    "Inclusive" kw governs whether to write lines which DO or DON'T meet the condition'''
    if out_txt_path is None:
        out_txt_path = in_txt_path.with_stem(f'{in_txt_path.stem}{"_" if postfix else ""}{postfix}')

    if (out_txt_path == in_txt_path):
        raise PermissionError(f'Attempting to overwrite {in_txt_path} with regex filter') # prevent write clash
    
    if (out_txt_path.suffix != in_txt_path.suffix):  # prevent file type conversion during transfer
        raise FileTypeError(f'Input and output file must have same extension (not {in_txt_path.suffix} and {out_txt_path.suffix})')

    with out_txt_path.open('w') as outfile: 
        with in_txt_path.open('r') as infile: # readfile is innermost in case error occurs during file read (caught by handler one level up)
            for line in infile:
                if (condition(line) == inclusive): # only write lines if (matching AND inclusive) OR (not matching AND exclusive)
                    outfile.write(line)

    if return_filtered_path:
        return out_txt_path