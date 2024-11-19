'''Tools for converting back and forth between integers and ordinal number words and prefices (https://en.wikipedia.org/wiki/Ordinal_numeral)'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

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