'''Base utilities and exceptions used throughout the bonding module''' # TODO : deprecate this module

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import logging
LOGGER = logging.getLogger(__name__)

from typing import Iterable, Union
from functools import reduce

from rdkit import Chem
from rdkit.Chem.rdchem import Mol, RWMol

from ..chemlabel import assign_contiguous_atom_map_nums


class BondOrderModificationError(Exception):
    '''Raised when an invalid RDKit bond modification is attempted'''
    pass

def combined_rdmol(*rdmols : Iterable[Mol], assign_map_nums : bool=True, editable : bool=True) -> Union[Mol, RWMol]:
    '''Merge any number of RDKit Mols into a single molecule with contiguous, non-overlapping atom map numbers'''
    if assign_map_nums: # assign contiguous, unique atom map numbers
        rdmols = assign_contiguous_atom_map_nums(*rdmols, in_place=False) 
    combo = reduce(Chem.CombineMols, rdmols) # combine into single Mol object to allow for bonding

    if editable:
        return Chem.RWMol(combo) # make combined Mol modifiable
    return combo