'''For mapping 1-to-1 between two allegedly identical molecules'''

from typing import Generator

from ..rdtypes import RDMol, RDAtom
from ..rderrors import SubstructMatchFailedError


# CUSTOM EXCEPTIONS
class MolSizeMismatchError(Exception):
    '''Raised when an operation is attempted on two molecules which were expected to have the same size but don't'''
    pass

# BIJECTIVE MATCHING METHODS
def bijective_atom_iter(rdmol_1 : RDMol, rdmol_2 : RDMol) -> Generator[tuple[RDAtom, RDAtom], None, None]:
    '''Takes two chemically identical molecules, matches corresponding atoms between them 1:1, and generates matching atom pairs
    Yields atoms in pairs in the same order as the molecules being matched were provided'''
    if rdmol_1.GetNumAtoms() != rdmol_2.GetNumAtoms():
        raise MolSizeMismatchError('Cannot generate bijection with Mols of two different sizes')
    
    atom_mapping = rdmol_1.GetSubstructMatch(rdmol_2) 
    if (not atom_mapping) or (len(atom_mapping) != rdmol_1.GetNumAtoms()): # TODO : verify whether the second condition necessary with the above size check
        raise SubstructMatchFailedError  # if no match or an incomplete match is found, raise Exception

    for rdatom_1_idx, rdatom_2 in zip(atom_mapping, rdmol_2.GetAtoms()):
        rdatom_1 = rdmol_1.GetAtomWithIdx(rdatom_1_idx)
        yield (rdatom_1, rdatom_2)