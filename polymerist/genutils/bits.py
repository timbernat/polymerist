'''For bitwise operations and conversions to/from bitstrings'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Union 


def int_to_bits(n : int, num_bits : int=8, clamp : bool=False, as_list : bool=False) -> Union[str, list[int]]:
    '''
    Convert an integer into a string of its bits, padded out to <num_bits> bits
    If clamp=True and the binary representation of <n> has more than <num_bits> bits, then leading bits will be discarded
    If as_list=True, will return as list of 0/1 ints; otherwise, will return as string
    '''
    bitstring = f'{n:0{num_bits}b}' # format into given number of bits with leading zeros - way easier than directly implementing bitshift-based operations
    if clamp:
        bitstring = bitstring[-num_bits:]

    if as_list:
        return [int(bit) for bit in bitstring] 
    return bitstring