'''For mapping 1-to-1 between two allegedly identical molecules'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Generator
from rdkit.Chem.rdchem import Atom, Mol

from ..rderrors import SubstructMatchFailedError


# CUSTOM EXCEPTIONS
class MolSizeMismatchError(Exception):
    '''Raised when an operation is attempted on two molecules which were expected to have the same size but don't'''
    pass

# BIJECTIVE MATCHING METHODS
def bijective_atom_id_iter(rdmol_1 : Mol, rdmol_2 : Mol) -> Generator[tuple[int, int], None, None]:
    '''Takes two chemically identical molecules, matches corresponding atoms between them 1:1, and generates matching atom id pairs
    Yields atoms in pairs in the same order as the molecules being matched were provided'''
    if rdmol_1.GetNumAtoms() != rdmol_2.GetNumAtoms():
        raise MolSizeMismatchError('Cannot generate bijection with Mols of two different sizes')
    
    atom_mapping = rdmol_1.GetSubstructMatch(rdmol_2) 
    if (not atom_mapping) or (len(atom_mapping) != rdmol_1.GetNumAtoms()): # TODO : verify whether the second condition necessary with the above size check
        raise SubstructMatchFailedError  # if no match or an incomplete match is found, raise Exception

    for rdatom_1_idx, rdatom_2 in zip(atom_mapping, rdmol_2.GetAtoms()):
        rdatom_2_idx = rdatom_2.GetIdx()
        yield (rdatom_1_idx, rdatom_2_idx)