'''For obtaining info from and for labelling individual RDKit Bonds'''

from typing import Any, Callable, Iterable
from itertools import combinations

from rdkit import Chem
from rdkit.Chem.rdchem import BondType

from ..rdtypes import RDMol, RDBond
from .molwise import atom_ids_by_map_nums
from ...genutils.typetools import Args, KWArgs
from ...genutils.iteration import sliding_window


# BOND ID QUERY FUNCTIONS
def get_bonded_pairs(rdmol : RDMol, *atom_ids : Iterable[int]) -> dict[int, tuple[int, int]]:
    '''Get bond and terminal atom indices of all bonds which exist between any pair of atoms in an indexed list'''
    res = {}
    
    atom_id_pairs = combinations(atom_ids, 2)
    for atom_id_pair in atom_id_pairs:
        bond = rdmol.GetBondBetweenAtoms(*atom_id_pair)
        if bond is not None:
            res[bond.GetIdx()] = atom_id_pair
    return res

def bond_ids_by_cond(rdmol : RDMol, bond_cond : Callable[[RDBond, Args, KWArgs], bool]) -> tuple[int]:
    '''Return IDs of all bonds which satisfy some binary condition'''
    return tuple(
        bond.GetIdx()
            for bond in rdmol.GetBonds()
                 if bond_cond(bond)
    )

def get_bonded_pairs_by_map_nums(rdmol : RDMol, *atom_map_nums : Iterable[int]) -> dict[int, tuple[int, int]]:
    '''Obtain bonded pair dict by atom map numbers instead of IDs'''
    return get_bonded_pairs(rdmol, *atom_ids_by_map_nums(rdmol, *atom_map_nums))

def get_bond_by_map_num_pair(rdmol : RDMol, map_num_pair : tuple[int, int]) -> RDBond:
    '''Get the index of a bond spanning a pair of atoms with given pair of atom map numbers'''
    return rdmol.GetBondBetweenAtoms(*atom_ids_by_map_nums(rdmol, *map_num_pair))

def get_bond_id_by_map_num_pair(rdmol : RDMol, map_num_pair : tuple[int, int]) -> RDBond:
    '''Get the index of a bond spanning a pair of atoms with given pair of atom map numbers'''
    return get_bond_by_map_num_pair(rdmol, map_num_pair).GetIdx()

def get_shortest_path_bonds(rdmol : RDMol, start_idx : int, end_idx : int) -> list[int]:
    '''Returns bond indices along shortest path between two atoms in a Mol'''
    return [
        rdmol.GetBondBetweenAtoms(*atom_id_pair).GetIdx()
            for atom_id_pair in sliding_window(Chem.GetShortestPath(rdmol, start_idx, end_idx), n=2)
    ]


