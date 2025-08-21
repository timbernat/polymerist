'''Utilities for validating, cleaning, and adding information into up SMILES and SMARTS strings'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Callable, TypeAlias, TypeVar
T = TypeVar('T')

from functools import wraps

from rdkit import Chem, RDLogger
from rdkit.Chem.rdmolops import SanitizeFlags, SanitizeMol, SANITIZE_ALL, SANITIZE_SETAROMATICITY


def suppress_rdkit_errors(func : Callable[..., T]) -> Callable[..., T]:
    '''Decorator to suppress RDKit error messages during function execution'''
    @wraps(func)
    def decorator(*args, **kwargs):
        RDLogger.DisableLog('rdApp.error')
        ret = func(*args, **kwargs)
        RDLogger.EnableLog('rdApp.error')
        
        return ret
    return decorator

# TYPING AND VALIDATION
Smiles : TypeAlias = str # purely for improving self-documentation of functions, no benefit to static type-checkers
Smarts : TypeAlias = str # purely for improving self-documentation of functions, no benefit to static type-checkers

@suppress_rdkit_errors
def is_valid_SMARTS(smarts : Smarts) -> bool:
    '''Check if SMARTS string is valid (according to RDKit)'''
    return (Chem.MolFromSmarts(smarts) is not None)

@suppress_rdkit_errors
def is_valid_SMILES(smiles : Smiles) -> bool:
    '''Check if SMARTS string is valid (according to RDKit)'''
    return (Chem.MolFromSmiles(smiles) is not None)

# CUSTOM EXCEPTIONS
class InvalidChemicalLineNotation(ValueError):
    '''Exception raised when a malformed chemical notation string is passed somewhere'''
    ...
    
## DEVNOTE: there are certainly more line notations out there; I'm just covering those actually used in the codebase here
class InvalidSMILES(InvalidChemicalLineNotation):
    '''Exception raised when a malformed SMILES string is passed somewhere'''
    ...

class InvalidSMARTS(InvalidChemicalLineNotation):
    '''Exception raised when a malformed SMARTS string is passed somewhere'''
    ...

class InvalidInChI(InvalidChemicalLineNotation):
    '''Exception raised when a malformed InChI string is passed somewhere'''
    ...
    

# CANONICALIZATION AND STRUCTURE EXPANSION
def canonical_SMILES_from_mol(mol : Chem.Mol) -> str:
    '''
    Cast Mol to a "canonical" SMILES format
    Mols with identical chemical structure should produce identical strings
    '''
    return Chem.CanonSmiles(Chem.MolToSmiles(mol, canonical=True))

def expanded_SMILES(
        smiles : str,
        assign_map_nums : bool=True,
        start_from : int=1,
        kekulize : bool=True,
        canonicalize : bool=True, # DEV: set to match legacy behavior
    ) -> str:
    '''
    Expands and clarifies the chemical information contained within a passed SMILES string
    namely explicit hydrogens and bond orders, and (optionally) kekulized aromatic bonds and atom map numbers
    '''
    if not is_valid_SMILES(smiles):
        raise InvalidSMILES(f'Passed string "{smiles}" cannot be interpreted as a valid SMILES pattern')
    
    rdmol = Chem.MolFromSmiles(smiles, sanitize=False)
    rdmol.UpdatePropertyCache() # inject valence and ring info without mangling from sanitization
    rdmol = Chem.AddHs(rdmol, addCoords=False)
    
    if assign_map_nums:
        for map_num, atom in enumerate(rdmol.GetAtoms(), start=start_from): # NOTE: deliberately did not use anything from rdutils.chemlabel here to avoid coupling
            atom.SetAtomMapNum(map_num) # NOTE that starting from anything below 1 will cause an atom somewhere to be mapped to 0 (i.e. not mapped)
    
    if kekulize:
        Chem.Kekulize(rdmol, clearAromaticFlags=True)
    Chem.SanitizeMol(rdmol)

    return Chem.MolToSmiles(rdmol, kekuleSmiles=kekulize, allBondsExplicit=True, allHsExplicit=True, canonical=canonicalize)