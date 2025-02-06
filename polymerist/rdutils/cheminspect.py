'''For inspecting what chemical information is present (or absent) in RDKit Mols'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from rdkit import Chem


def all_Hs_are_explicit(mol : Chem.Mol) -> bool:
    '''Determine whether an RDKit Mol has all hydrogens explicitly present'''
    return all(
        atom.GetNumImplicitHs() == 0
            for atom in mol.GetAtoms()
    )
    
def has_aromatic_bonds(mol : Chem.Mol) -> bool:
    '''Determine whether an RDKit Mol has any aromatic bonds present'''
    return any(
        bond.GetBondType() == Chem.BondType.AROMATIC
            for bond in mol.GetBonds()
    )