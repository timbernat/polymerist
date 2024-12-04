'''For searching and replacing through strings and text files'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'


def uniquify_str(string : str, preserve_order : bool=True) -> str:
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
            uniquify_str("balaclava", preserve_order=False) -> "bcavl"
            uniquify_str("balaclava", preserve_order=True) -> "balcv"
        
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