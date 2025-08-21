'''Encoding, hashing, and conversion of string to and from various formats'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import re
import hashlib, base64


# HASHING 
def hash_as_alphanumeric(string : str, hash_algorithm : str='md5', str_encoding : str='utf-8', strip_b64_padding : bool=True) -> str:
    '''Map string to hash text which contains only alphanumeric + dash characters'''
    if (hash_algorithm not in hashlib.algorithms_guaranteed):
        raise KeyError(f'Invalid hash algorithm "{hash_algorithm}". Supported algorithms are: {", ".join(hashlib.algorithms_guaranteed)}')
    hash_funct = getattr(hashlib, hash_algorithm)

    hashbytes = hash_funct(string.encode(str_encoding)).digest()
    hashtext  = base64.urlsafe_b64encode(hashbytes).decode(str_encoding)
    if strip_b64_padding:
        hashtext = hashtext.rstrip('=') # remove padding to avoid possibly invalid characters (and ugliness)

    return hashtext
hash_as_alphanum = hash_as_alphanumeric # alias for convenience

# TYPE COERCION
INT_REGEX = re.compile(
    r'''
    (?!.*_$)        # lookahead and fail is line ends contains any number of characters ending with underscore...
    ^               # otherwise match line start
    [-+]?           # an optional plus or minus sign
    (               # then either: 
        [1-9]       # 1) any nonzero digit (specifically NOT 0 or an underscore)...
        [\d_]*      # ...followed by any number (including 0) of digits or underscores
        |0          # 2) 0 by itself as a special case
    )$              # then finally the end of line and nothing more
    ''',
    flags=re.VERBOSE,
)
def representable_as_int(string : str) -> bool:
    '''Check if a string corresponds to (i.e. can be represented as) a well-defined Python int'''
    # try:
    #     _ = int(string) # overly inclusive to floats etc. which can be cast as ints, even incases where we want to enforce "int-itude"
    # except ValueError:
    #     return False
    # else:
    #     return True
    return re.match(INT_REGEX, string) is not None