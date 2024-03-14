'''SMILES and SMARTS primitives and functions for validation'''

from rdkit import Chem
from rdkit.Chem.rdchem import BondType


BOND_SMILES_BY_ORDER = { # concrete bond objects (since these can't be directly instantiated from bond order in Python)
    BondType.SINGLE    : Chem.BondFromSmiles('-'),
    BondType.DOUBLE    : Chem.BondFromSmiles('='),
    BondType.TRIPLE    : Chem.BondFromSmiles('#'),
    BondType.QUADRUPLE : Chem.BondFromSmiles('$'),
    BondType.AROMATIC  : Chem.BondFromSmiles(':'),
}

BOND_SMARTS_BY_ORDER = { # concrete bond objects (since these can't be directly instantiated from bond order in Python)
    BondType.SINGLE    : Chem.BondFromSmarts('-'),
    BondType.DOUBLE    : Chem.BondFromSmarts('='),
    BondType.TRIPLE    : Chem.BondFromSmarts('#'),
    BondType.QUADRUPLE : Chem.BondFromSmarts('$'),
    BondType.AROMATIC  : Chem.BondFromSmarts(':'),
}

def is_valid_SMARTS(smarts : str) -> bool:
    '''Check if SMARTS string is valid (according to RDKit)'''
    return (Chem.MolFromSmarts(smarts) is not None)

def is_valid_SMILES(smiles : str) -> bool:
    '''Check if SMARTS string is valid (according to RDKit)'''
    return (Chem.MolFromSmiles(smiles) is not None)