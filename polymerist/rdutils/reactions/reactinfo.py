'''Classes for encapsulating conserved info about molecular components during a reaction'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Union, Optional, Iterable, Sequence
from dataclasses import dataclass
from enum import StrEnum, auto

from rdkit.Chem import Mol, BondType

from ..rdprops import copy_rdobj_props


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
    
    @staticmethod
    def apply_atom_info_to_product(
            product : Mol,
            product_atom_infos : Iterable['AtomTraceInfo'],
            reactants : Sequence[Mol],
            apply_map_labels : bool=True,
        ) -> None:
        '''Transfer props and (if requested) map number information from atoms in reactant Mols to their corresponding atoms in a product Mol
        Acts in-place on the "product" Mol instance'''
        for atom_info in product_atom_infos:
            product_atom = product.GetAtomWithIdx(atom_info.product_atom_idx)
            assert product_atom.HasProp('old_mapno') # precisely the mapped atoms in the reaction template will have this property set on the product 
            assert product_atom.GetIntProp('old_mapno') == atom_info.map_number # sanity check to make sure atom map identity has been preserved

            corresp_reactant_atom = reactants[atom_info.reactant_idx].GetAtomWithIdx(atom_info.reactant_atom_idx) # product_atom.GetIntProp('react_atom_idx'))
            copy_rdobj_props(from_rdobj=corresp_reactant_atom, to_rdobj=product_atom) # clone props from corresponding atom in reactant
            
            if apply_map_labels:
                product_atom.SetAtomMapNum(atom_info.map_number)
                product_atom.ClearProp('old_mapno')

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
    
    @staticmethod
    def apply_bond_info_to_product(
            product : Mol,
            product_bond_infos : Iterable['BondTraceInfo'],
        ) -> None:
        '''Mark any changed bonds with bond props and clean up bond type info in places where bonds get modified
        Acts in-place on the "product" Mol instance'''
        for bond_info in product_bond_infos:
            product_bond = product.GetBondWithIdx(bond_info.product_bond_idx)
            
            if bond_info.bond_change_type == BondChange.ADDED: # explicictly labl new or changed bonds
                product_bond.SetProp(BOND_CHANGE_PROPNAME, BondChange.ADDED)
    
            if bond_info.bond_change_type == BondChange.MODIFIED:
                product_bond.SetProp(BOND_CHANGE_PROPNAME, BondChange.MODIFIED)
                assert(product_bond.GetBeginAtom().HasProp('_ReactionDegreeChanged')) # double check reaction agrees the bond has changed
                assert(product_bond.GetEndAtom().HasProp('_ReactionDegreeChanged')) 

                if bond_info.final_bond_type not in (BondType.ZERO, BondType.UNSPECIFIED): # if an explicit bond type is defined in the template...
                    product_bond.SetBondType(bond_info.final_bond_type) # ...set the product's bond type to what it *should* be from the reaction schema
