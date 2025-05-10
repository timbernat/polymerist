'''For assigning, transferring, and removing properties of RDKit Atoms'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Callable, Iterable, Optional, TypeVar
T = TypeVar('T')

from rdkit.Chem import Atom, Mol

from ._rdprops import RDPROP_GETTERS
from ...genutils.decorators.functional import optional_in_place


# ATOM PROPERTY INSPECTION
def atom_ids_with_prop(rdmol : Mol, prop_name : str) -> list[int]:
    '''Returns list of indices of atoms which have a particular property assigned'''
    return [
        atom.GetIdx()
            for atom in rdmol.GetAtoms()
                if atom.HasProp(prop_name)
    ]

def aggregate_atom_prop(rdmol : Mol, prop : str, prop_type : T=str) -> dict[int, T]:
    '''Collects the values of a given Prop across all atoms in an RDKit molecule'''
    getter_type = RDPROP_GETTERS[prop_type]
    return {
        atom_idx : getattr(rdmol.GetAtomWithIdx(atom_idx), getter_type)(prop) # invoke type-appropriate getter on atom, with the name of the desired property
            for atom_idx in atom_ids_with_prop(rdmol, prop_name=prop)
    }
    
@optional_in_place
def clear_atom_props(rdmol : Mol) -> None: # TODO: add support for clearing only specific named/private props?
    '''Wipe properties of all atoms in a molecule'''
    for atom in rdmol.GetAtoms():
        for prop_name in atom.GetPropNames():
            atom.ClearProp(prop_name)

# ATOM PROP ANNOTATION
@optional_in_place
def annotate_atom_prop(rdmol : Mol, prop : str, prop_type : T=str, annotate_precision : Optional[int]=None) -> None:
    '''Labels the desired Prop for all atoms in a Mol which have it'''
    getter_type = RDPROP_GETTERS[prop_type]
    for atom in rdmol.GetAtoms():
        prop_val = getattr(atom, getter_type)(prop) # invoke type-appropriate getter on atom, with the name of the desired property
        
        if hasattr(prop_val, '__round__') and annotate_precision is not None: # only round on roundable objects, and only when
            prop_val = round(prop_val, annotate_precision)
        atom.SetProp('atomNote', str(prop_val)) # need to convert to string, as double is susceptible to float round display errors (shows all decimal places regardless of rounding)
    
@optional_in_place
def annotate_atom_ids(rdmol : Mol, atom_id_remap : Optional[dict[int, int]]=None) -> None:
    '''Draws atom indices over their positions when displaying a Mol. 
    Can optionally provide a dict mapping atom indices to some other integers'''
    if atom_id_remap is None:
        atom_id_remap = {} # avoid mutable default

    for atom in rdmol.GetAtoms():
        atom.SetIntProp('atomNote', atom_id_remap.get(atom.GetIdx(), atom.GetIdx())) # check if map value exists, if not default to index

@optional_in_place
def clear_atom_annotations(rdmol : Mol) -> None:
    '''Removes atom annotations over their positions when displaying a Mol'''
    for atom in rdmol.GetAtoms():
        atom.ClearProp('atomNote')
        
@optional_in_place
def label_linkers(rdmol : Mol, label_props : Iterable[str]=None, naming_funct : Optional[Callable[[int], str]]=None) -> None:
    '''Labels wild-type ("*") atoms in a Mol for display'''
    if label_props is None:
        label_props = ('_displayLabel',) # by default, only set labels for display
    
    if naming_funct is None:
        naming_funct = lambda map_num : f'R<sub>{map_num or ""}<sub>' # will be just "R" if atom is unmapped
    
    for atom in rdmol.GetAtoms():
        if atom.GetAtomicNum() != 0:
            continue # skip explicit atoms
        
        map_num = atom.GetAtomMapNum()
        for label_prop in label_props:
            atom.SetProp(label_prop, naming_funct(map_num))
label_R_groups = label_linkers
