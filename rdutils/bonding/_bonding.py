'''Base utilities and exceptions used throughout the bonding module'''

import logging
LOGGER = logging.getLogger(__name__)

from typing import Iterable, Union
from functools import reduce
from rdkit import Chem

from ..rdtypes import RWMol, RDMol
from ..labeling import molwise


# BOND REFERENCE
class BondOrderModificationError(Exception):
    '''Raised when an invalid RDBond modification is attempted'''
    pass

def combined_rdmol(*rdmols : Iterable[RDMol], assign_map_nums : bool=True, editable : bool=True) -> Union[RDMol, RWMol]:
    '''Merge any number of RDMols into a single molecule with contiguous, non-overlapping atom map numbers'''
    if assign_map_nums: # assign contiguous, unique atom map numbers
        rdmols = molwise.assign_contiguous_atom_map_nums(*rdmols, in_place=False) 
    combo = reduce(Chem.CombineMols, rdmols) # combine into single Mol object to allow for bonding

    if editable:
        return Chem.RWMol(combo) # make combined Mol modifiable
    return combo