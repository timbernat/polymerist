'''For searching and replacing through strings and text files'''

from typing import Callable, Optional
from pathlib import Path

from ..fileutils.extensions import FileTypeError


def shortest_repeating_substring(string : str) -> str:
    '''Return the shortest substring such that the passed string can be written as some number of repeats (including 1) of the substring
    Will return the original string if no simpler decomposition exists'''
    i = (2*string).find(string, 1, -1) # check if string matches itself in a cycle in non-trivial way (i.e more than just the two repeats)
    return string if (i == -1) else string[:i]

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