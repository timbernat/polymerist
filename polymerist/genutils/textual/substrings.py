'''For identifying and concatenating substrings of other strings with unique properties'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'


def unique_string(string : str, preserve_order : bool=True) -> str:
    '''
    Accepts a string and returns another string containing
    only the UNIQUE characters in the origin string
    
    Can specify whether order is important with the "preserve_order" keyword
    
    Parameters
    ----------
    string : str
        An arbitrary string on wants the unique characters from
    preserve_order : bool, default True
        Whether or not to keep the unique characters in the order they are found
        For example: 
            unique_string("balaclava", preserve_order=False) -> "bcavl"
            unique_string("balaclava", preserve_order=True) -> "balcv"
        
    Returns
    -------
    uniquified_str : str
        Another string containing only the unique characters in "string"
        Order depends on the value of the "preserve_order" parameter
    '''
    if not preserve_order:
        unique_chars = set(string)
    else:
        unique_chars = []
        for char in string:
            if char not in unique_chars:
                unique_chars.append(char)
    
    return ''.join(unique_chars)

def shortest_repeating_substring(string : str) -> str:
    '''Return the shortest substring such that the passed string can be written as some number of repeats (including 1) of the substring
    Will return the original string if no simpler decomposition exists'''
    i = (2*string).find(string, 1, -1) # check if string matches itself in a cycle in non-trivial way (i.e more than just the two repeats)
    return string if (i == -1) else string[:i]

def repeat_string_to_length(string : str, target_length : int, joiner : str='') -> str:
    '''
    Takes a string and repeats it cyclically to produce another string of a given length
    The number of times the original string occurs in the new string may be fractional
    for example:
    >> repeat_string_to_length("CAT", 6) -> "CATCAT"
    >> repeat_string_to_length("BACA", 10) -> "BACABACABA"
    
    Parameters
    ----------
    string : str
        An arbitrary string to repeat
    target_length : int
        The length of the final desired string
        This does NOT have to be an integer multiple of the length of "string"
            E.g. repeat_string_to_length("BACA", 10) -> "BACABACABA"
        Nor does it have to be greater than the length of "string"
            E.g. repeat_string_to_length("BACA", 3) -> "BAC"
            
    Returns
    -------
    rep_string : str
        A new string which has the desired target length and consists of cycles of the initial string
    '''
    if not string:
        raise ValueError(f'Cannot generate nonempty string from any amount of repeats of the empty string')
    if not isinstance(target_length, int):
        raise TypeError(f'Only integer target string lengths are allowed, not non-integer type "{type(target_length).__name__}"')
    if target_length < 0:
        raise IndexError(f'Cannot generate a string of negative length (requested length of {target_length} character(s))')
    
    num_str_reps, num_extra_chars = divmod(target_length, len(string))
    remainder = (string[:num_extra_chars],) if num_extra_chars else () # empty container avoids extra joiner at end when remainder string is empty
    
    return joiner.join(num_str_reps*(string,) + remainder) # tuples here are ~2 OOM faster than moral equivalent with lists
    