'''For mapping 1-to-1 between pairs of identical molecules'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Generator
from rdkit.Chem.rdchem import Atom, Mol


# BIJECTIVE ATOM MAPPING
class SubstructMatchFailedError(Exception):
    '''Raised when molecule graph isomorphism match does not form a cover'''
    pass

class MolSizeMismatchError(Exception):
    '''Raised when an operation is attempted on two molecules which were expected to have the same size but don't'''
    pass

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

# ATOM-BY-ATOM PROPERTY COMPARISON BETWEEN MOLECULES        
def difference_rdmol(rdmol_1 : Mol, rdmol_2 : Mol, prop : str='PartialCharge', remove_map_nums : bool=True) -> Mol:
    '''
    Takes two RDKit Mols (presumed to have the same structure and atom map numbers) and the name of a property 
    whose partial charges are the differences betwwen the two Mols' charges (atomwise)
    
    Assumes that the property in question is numeric (i.e. can be interpreted as a float)
    '''
    diff_mol = Mol(rdmol_1) # duplicate first molecule as template
    all_deltas = []
    for rdatom_1_idx, rdatom_2_idx in bijective_atom_id_iter(rdmol_1, rdmol_2):
        rdatom_1 = rdmol_1.GetAtomWithIdx(rdatom_1_idx)
        rdatom_2 = rdmol_2.GetAtomWithIdx(rdatom_2_idx)
        diff_atom = diff_mol.GetAtomWithIdx(rdatom_1_idx) # same index, since it is a deep copy

        delta = rdatom_1.GetDoubleProp(prop) - rdatom_2.GetDoubleProp(prop)
        diff_atom.SetDoubleProp(f'Delta{prop}', delta)
        all_deltas.append(delta)

        diff_atom.ClearProp(prop) # reset property value from original atom copy to avoid confusion
        if remove_map_nums:
            diff_atom.ClearProp('molAtomMapNumber') # Remove atom map num for greater visual clarity when drawing

    diff_mol.ClearProp(prop) # reset property value from original mol copy to avoid confusion
    diff_mol.SetProp(f'Delta{prop}s', str(all_deltas)) # label stringified version of property list (can be de-stringified via ast.literal_eval)
    diff_mol.SetDoubleProp(f'Delta{prop}Min', min(all_deltas)) # label minimal property value for ease of reference
    diff_mol.SetDoubleProp(f'Delta{prop}Max', max(all_deltas)) # label maximal property value for ease of reference

    return diff_mol