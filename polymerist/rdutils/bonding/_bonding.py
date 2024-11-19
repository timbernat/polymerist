'''Base utilities and exceptions used throughout the bonding module''' # TODO : deprecate this module

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import logging
LOGGER = logging.getLogger(__name__)

from typing import Iterable, Union
from functools import reduce

from rdkit import Chem
from rdkit.Chem.rdchem import Mol, RWMol

from ..labeling import molwise


# BOND REFERENCE
def combined_rdmol(*rdmols : Iterable[Mol], assign_map_nums : bool=True, editable : bool=True) -> Union[Mol, RWMol]:
    '''Merge any number of RDKit Mols into a single molecule with contiguous, non-overlapping atom map numbers'''
    if assign_map_nums: # assign contiguous, unique atom map numbers
        rdmols = molwise.assign_contiguous_atom_map_nums(*rdmols, in_place=False) 
    combo = reduce(Chem.CombineMols, rdmols) # combine into single Mol object to allow for bonding

    if editable:
        return Chem.RWMol(combo) # make combined Mol modifiable
    return combo