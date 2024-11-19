'''Encoding, hashing, and conversion of string to and from various formats'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import hashlib, base64


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