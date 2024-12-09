'''SMILES and SMARTS primitives and functions for validation'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import TypeAlias

from rdkit import Chem
from rdkit.Chem.rdchem import BondType


# VALIDATION
Smiles : TypeAlias = str # purely for improving self-documentation of functions, no benefit to static type-checkers
Smarts : TypeAlias = str # purely for improving self-documentation of functions, no benefit to static type-checkers

def is_valid_SMARTS(smarts : Smarts) -> bool:
    '''Check if SMARTS string is valid (according to RDKit)'''
    return (Chem.MolFromSmarts(smarts) is not None)

def is_valid_SMILES(smiles : Smiles) -> bool:
    '''Check if SMARTS string is valid (according to RDKit)'''
    return (Chem.MolFromSmiles(smiles) is not None)

# BOND PRIMITIVES AND RELATED OBJECTS
BOND_PRIMITIVES = '~-=#$:'
BOND_PRIMITIVES_FOR_REGEX = r'[~\-=#$:]' # any of the SMARTS bond primitive chars, with a space to differentiate single-bond hyphen for the regex range char
BOND_INITIALIZERS = {
    'SMILES' : (Chem.Bond     , Chem.BondFromSmiles),
    'SMARTS' : (Chem.QueryBond, Chem.BondFromSmarts),
}
    
for in_line_fmt, (rd_type, rd_initializer) in BOND_INITIALIZERS.items():
    in_line_fmt = in_line_fmt.upper()
    rd_prefix   = rd_type.__name__.upper()

    globals()[f'BONDTYPE_BY_BOND_{in_line_fmt}'] = type_by_bonds  = {}
    globals()[f'ORDER_BY_BOND_{in_line_fmt}'   ] = order_by_bonds = {}
    globals()[f'BOND_{in_line_fmt}_BY_BONDTYPE'] = bonds_by_type  = {}
    globals()[f'BOND_{in_line_fmt}_BY_ORDER'   ] = bonds_by_order = {}
    globals()[f'RDKIT_{rd_prefix}S_BY_BONDTYPE'] = rdbonds_by_type  = {}
    globals()[f'RDKIT_{rd_prefix}S_BY_ORDER'   ] = rdbonds_by_order = {}
    
    for prim_str in BOND_PRIMITIVES:
        rd_bond = rd_initializer(prim_str)
        if (rd_bond is not None) and (type(rd_bond) == rd_type):
            bondtype, order = rd_bond.GetBondType(), rd_bond.GetBondTypeAsDouble()
            type_by_bonds[prim_str]   = bondtype
            order_by_bonds[prim_str]  = order
            bonds_by_type[bondtype]   = prim_str
            bonds_by_order[order]     = prim_str
            rdbonds_by_type[bondtype] = rd_bond
            rdbonds_by_order[order]   = rd_bond
