'''Utilities for validating, cleaning, and adding information into up SMILES and SMARTS strings'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import TypeAlias

from rdkit import Chem
from rdkit.Chem.rdmolops import SanitizeFlags, SanitizeMol, SANITIZE_ALL, SANITIZE_SETAROMATICITY

from ..rdutils.labeling.molwise import assign_ordered_atom_map_nums


# TYPING AND VALIDATION
Smiles : TypeAlias = str # purely for improving self-documentation of functions, no benefit to static type-checkers
Smarts : TypeAlias = str # purely for improving self-documentation of functions, no benefit to static type-checkers

def is_valid_SMARTS(smarts : Smarts) -> bool:
    '''Check if SMARTS string is valid (according to RDKit)'''
    return (Chem.MolFromSmarts(smarts) is not None)

def is_valid_SMILES(smiles : Smiles) -> bool:
    '''Check if SMARTS string is valid (according to RDKit)'''
    return (Chem.MolFromSmiles(smiles) is not None)

# SANITIZATION (IN THE RDKit SENSE)
# TODO: add decorator for sanitization flag inject of functions which return mols
SANITIZE_AS_KEKULE = (Chem.SANITIZE_ALL & ~Chem.SANITIZE_SETAROMATICITY) # sanitize everything EXCEPT reassignment of aromaticity

# CANONICALIZATION
def canonicalize_mol(mol : Chem.Mol) -> str:
    '''
    Cast Mol to a "canonical" SMILES format -
    Mols with identical chemical structure should produce identical strings
    '''
    return Chem.CanonSmiles(Chem.MolToSmiles(mol, canonical=True))

# STRUCTURE EXPANSION
def expanded_SMILES(
        smiles : str,
        assign_map_nums : bool=True,
        start_from : int=1,
        kekulize : bool=True,
    ) -> str:
    '''
    Expands and clarifies the chemical information contained within a passed SMILES string
    namely explicit hydrogens and bond orders, and (optionally) kekulized aromatic bonds and atom map numbers
    '''
    assert(is_valid_SMILES(smiles))
    
    rdmol = Chem.MolFromSmiles(smiles, sanitize=True)
    rdmol = Chem.AddHs(rdmol, addCoords=True)
    if assign_map_nums:
        rdmol = assign_ordered_atom_map_nums(rdmol, start_from=start_from)
    
    if kekulize:
        Chem.Kekulize(rdmol, clearAromaticFlags=True)
    Chem.SanitizeMol(rdmol)

    return Chem.MolToSmiles(rdmol, kekuleSmiles=kekulize, allBondsExplicit=True, allHsExplicit=True)