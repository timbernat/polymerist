'''Classes for encapsulating conserved info about molecular components during a reaction'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Union, Optional, Iterable, Sequence
from dataclasses import dataclass
from enum import StrEnum, auto

from rdkit.Chem import Mol, BondType


REACTANT_INDEX_PROPNAME : str = 'reactant_idx' # name of the atom property to assign reactant template indices to
BOND_CHANGE_PROPNAME : str = 'bond_change'  # name of bond property to set on bonds to indicate they have changed in a reaction

class BondChange(StrEnum):
    '''For indicating how a bond which changed in a reaction was altered'''
    ADDED = auto()
    DELETED = auto()
    MODIFIED = auto() # specifically, when bond order is modified but the bond persists
    UNCHANGED = auto()
    
@dataclass(frozen=True)
class AtomTraceInfo:
    '''For encapsulating information about the origin and destination of a mapped atom, traced through a reaction'''
    map_number : int
    reactant_idx      : int # index of the reactant template within a reaction in which the atom occurs
    reactant_atom_idx : int # index of the target atom WITHIN the above reactant template
    product_idx       : int # index of the product template within a reaction in which the atom occurs
    product_atom_idx  : int # index of the target atom WITHIN the above product template

@dataclass(frozen=True)
class BondTraceInfo:
    '''For encapsulating information about bonds which are between mapped atoms and which change during a reaction'''
    map_nums : Union[tuple[int, int], frozenset[int]] # map numbers of the pair of atoms the bond connects
    # NOTE: reactant index doesn't make much sense, since the atoms the bond spans might have comes from two distinct reactant templates
    # product index is also debatable, since deleted bonds may place atoms into separate products in general
    product_idx       : Optional[int] # index of the reactant template within a reaction in which the modified bond occurs
    product_bond_idx  : Optional[int] # index of the target bond WITHIN the above product template
    bond_change_type  : Union[str, BondChange]
    initial_bond_type : Union[float, BondType] # bond order in the reactant (i.e. BEFORE the change)
    final_bond_type   : Union[float, BondType] # bond order in the product (i.e. AFTER the change)
