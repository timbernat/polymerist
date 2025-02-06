'''For obtaining info from and for labelling individual RDKit Bonds'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Iterable, Optional, ParamSpec, Union
Params = ParamSpec('Params')

from itertools import combinations
from rdkit.Chem.rdchem import Bond, BondType, Mol

from .molwise import get_atom_ids_by_map_nums


# BOND ID QUERYING    
def get_bonded_pairs(rdmol : Mol, *atom_ids : Iterable[int]) -> dict[int, tuple[int, int]]:
    '''Get bond and (begin,end) atom indices of all bonds which exist between any pair of atoms in an indexed list'''
    return {
        bond.GetIdx() : atom_id_pair
            for atom_id_pair in combinations(atom_ids, 2)
                if (bond := rdmol.GetBondBetweenAtoms(*atom_id_pair)) is not None
    }

def get_bonded_pairs_by_map_nums(rdmol : Mol, *atom_map_nums : Iterable[int]) -> dict[int, tuple[int, int]]:
    '''Obtain bonded pair dict by atom map numbers instead of IDs'''
    return get_bonded_pairs(rdmol, *get_atom_ids_by_map_nums(rdmol, *atom_map_nums))

def get_bond_by_map_num_pair(rdmol : Mol, map_num_pair : tuple[int, int], as_bond : bool=True) -> Optional[Union[int, Bond]]:
    '''
    Get the bond spanning a pair of atoms with given pair of atom map numbers
    Returns the RDkit.Bond object if as_bond=True, and the index of the bond if as_bond=False
    
    If no bond exists between the atoms, will return None regardless of the value of "as_bond"
    '''
    bond = rdmol.GetBondBetweenAtoms(*get_atom_ids_by_map_nums(rdmol, *map_num_pair))
    if (not as_bond) and (bond is not None):
        return bond.GetIdx()
    return bond # returns bond or, implicitly, NoneType if no bond is found
