'''For assigning, transferring, and removing properties of RDKit Atoms'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Callable, Generator, Optional, TypeVar
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
def clear_atom_props(rdmol : Mol) -> None:
    '''Wipe properties of all atoms in a molecule'''
    for atom in rdmol.GetAtoms():
        for prop_name in atom.GetPropNames():
            atom.ClearProp(prop_name)

# ATOM PROP ANNOTATION
@optional_in_place
def annotate_atom_ids(rdmol : Mol, atom_id_remap : Optional[dict[int, int]]=None) -> None:
    '''Draws atom indices over their positions when displaying a Mol. 
    Can optionally provide a dict mapping atom indices to some other integers'''
    if atom_id_remap is None:
        atom_id_remap = {} # avoid mutable default

    for atom in rdmol.GetAtoms():
        atom.SetIntProp('atomNote', atom_id_remap.get(atom.GetIdx(), atom.GetIdx())) # check if map value exists, if not default to index

            
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
def clear_atom_annotations(rdmol : Mol) -> None:
    '''Removes atom annotations over their positions when displaying a Mol'''
    for atom in rdmol.GetAtoms():
        atom.ClearProp('atomNote')

# ATOM NEIGHBOR SEARCH
def _get_atom_neighbors_by_condition_factory(condition : Callable[[Atom], bool]) -> Callable[[Atom], Generator[Atom, None, None]]:
    '''Factory function for generating neighbor-search functions over Atoms by a boolean condition'''
    def neighbors_by_condition(atom : Atom) -> Generator[Atom, None, None]:
        '''Generate all neighboring atoms satisfying a condition'''
        for nb_atom in atom.GetNeighbors():
            if condition(nb_atom):
                yield nb_atom
    return neighbors_by_condition

def _has_atom_neighbors_by_condition_factory(condition : Callable[[Atom], bool]) -> Callable[[Atom], bool]:
    '''Factory function for generating neighbor-search functions over Atoms by a boolean condition'''
    def has_neighbors_by_condition(atom : Atom) -> bool:
        '''Identify if any neighbors of an atom satisfy some condition'''
        return any(
            condition(nb_atom)
                for nb_atom in atom.GetNeighbors()
        )
    return has_neighbors_by_condition