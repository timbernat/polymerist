'''Simplified typehints for RDKit-related functionality'''

from typing import Any, TypeAlias, Union
from rdkit import Chem

from ..genutils.typetools.categorical import _union_member_factory

# DEFINING RDKit TYPES
RDMol  = Chem.rdchem.Mol
RWMol  = Chem.rdchem.RWMol
RDAtom = Chem.rdchem.Atom
RDBond = Chem.rdchem.Bond

RDObj : TypeAlias = Union[RDMol, RWMol, RDAtom, RDBond]
isrdobj = _union_member_factory(RDObj, 'RDObj')