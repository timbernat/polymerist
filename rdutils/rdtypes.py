'''Simplified typehints for RDKit-related functionality'''

from typing import Any, TypeAlias, Union
from rdkit import Chem


# DEFINING RDKit TYPES
RDMol  = Chem.rdchem.Mol
RWMol  = Chem.rdchem.RWMol
RDAtom = Chem.rdchem.Atom
RDBond = Chem.rdchem.Bond

RDObj : TypeAlias = Union(RDMol, RWMol, RDAtom, RDBond)

# TYPECHECKING
def isrdobj(var : Any) -> bool:
	'''Check if a variable is and RDKit object'''
	return isinstance(var, RDObj.__args__)