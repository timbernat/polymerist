'''For assigning, transferring, and removing properties of RDKit Bonds'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Optional, TypeVar
T = TypeVar('T')

from rdkit.Chem import Bond, Mol

from ._rdprops import RDPROP_GETTERS
from ...genutils.decorators.functional import optional_in_place


# BOND PROPERTY INSPECTION
def bond_ids_with_prop(rdmol : Mol, prop_name : str) -> list[int]:
    '''Returns list of indices of bonds which have a particular property assigned'''
    return [
        bond.GetIdx()
            for bond in rdmol.GetBonds()
                if bond.HasProp(prop_name)
    ]

def aggregate_bond_prop(rdmol : Mol, prop : str, prop_type : T=str) -> dict[int, T]:
    '''Collects the values of a given Prop across all bonds in an RDKit molecule'''
    getter_type = RDPROP_GETTERS[prop_type]
    return {
        bond_idx : getattr(rdmol.GetBondWithIdx(bond_idx), getter_type)(prop) # invoke type-appropriate getter on bond, with the name of the desired property
            for bond_idx in bond_ids_with_prop(rdmol, prop_name=prop)
    }
    
@optional_in_place
def clear_bond_props(rdmol : Mol) -> None:
    '''Wipe properties of all bonds in a molecule'''
    for bond in rdmol.GetBonds():
        for prop_name in bond.GetPropNames():
            bond.ClearProp(prop_name)
            
# BOND PROP ANNOTATION
@optional_in_place
def annotate_bond_prop(rdmol : Mol, prop : str, prop_type : T=str, annotate_precision : Optional[int]=None) -> None:
    '''Labels the desired Prop for all bonds in a Mol which have it'''
    getter_type = RDPROP_GETTERS[prop_type]
    for bond in rdmol.GetBonds():
        prop_val = getattr(bond, getter_type)(prop) # invoke type-appropriate getter on bond, with the name of the desired property
        
        if hasattr(prop_val, '__round__') and annotate_precision is not None: # only round on roundable objects, and only when
            prop_val = round(prop_val, annotate_precision)
        bond.SetProp('bondNote', str(prop_val)) # need to convert to string, as double is susceptible to float round display errors (shows all decimal places regardless of rounding)
    
@optional_in_place
def annotate_bond_ids(rdmol : Mol, bond_id_remap : Optional[dict[int, int]]=None) -> None:
    '''Draws bond indices over their positions when displaying a Mol. 
    Can optionally provide a dict mapping bond indices to some other integers'''
    if bond_id_remap is None:
        bond_id_remap = {} # avoid mutable default

    for bond in rdmol.GetBonds():
        bond.SetIntProp('bondNote', bond_id_remap.get(bond.GetIdx(), bond.GetIdx())) # check if map value exists, if not default to index

@optional_in_place
def clear_bond_annotations(rdmol : Mol) -> None:
    '''Removes bond annotations over their positions when displaying a Mol'''
    for bond in rdmol.GetBonds():
        bond.ClearProp('bondNote')
        