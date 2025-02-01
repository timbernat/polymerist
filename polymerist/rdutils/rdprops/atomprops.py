'''For assigning, transferring, and removing properties of RDKit Atoms'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Callable, Generator, Optional, TypeVar
T = TypeVar('T')   # generic type for AtomProp attributes

from rdkit.Chem import Atom, Mol

from .rdprops import RDPROP_GETTERS
from ...genutils.decorators.functional import optional_in_place


# REFERENCE FOR "MAGIC" ATOM PROP KEYS (https://www.rdkit.org/docs/RDKit_Book.html#atom)
RDATOM_MAGIC_PROPS = {
    '_CIPCode'               : 'the CIP code (R or S) of the atom',
    '_CIPRank'               : 'the integer CIP rank of the atom',
    '_ChiralityPossible'     : 'set if an atom is a possible chiral center',
    '_MolFileRLabel'         : 'integer R group label for an atom, read from/written to CTABs.',
    '_ReactionDegreeChanged' : 'set on an atom in a product template of a reaction if its degree changes in the reaction',
    '_protected'             : 'atoms with this property set will not be considered as matching reactant queries in reactions',
    'dummyLabel'             : '(on dummy atoms) read from/written to CTABs as the atom symbol',
    'molAtomMapNumber'       : 'the atom map number for an atom, read from/written to SMILES and CTABs',
    'molfileAlias'           : 'the mol file alias for an atom (follows A tags), read from/written to CTABs',
    'molFileValue'           : 'the mol file value for an atom (follows V tags), read from/written to CTABs',
    'molFileInversionFlag'   : 'used to flag whether stereochemistry at an atom changes in a reaction, read from/written to CTABs, determined automatically from SMILES',
    'molRxnComponent'        : 'which component of a reaction an atom belongs to, read from/written to CTABs',
    'molRxnRole'             : 'which role an atom plays in a reaction (1=Reactant, 2=Product, 3=Agent), read from/written to CTABs',
    'smilesSymbol'           : 'determines the symbol that will be written to a SMILES for the atom',
}

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

# ATOM PROPERTY INSPECTION
def atom_ids_with_prop(rdmol : Mol, prop_name : str) -> list[int]:
    '''Returns list of atom IDs of atom which have a particular property assigned'''
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
def annotate_atom_prop(rdmol : Mol, prop : str, prop_type : T=str, annotate_precision : Optional[int]=None) -> None:
    '''Labels the desired Prop for all atoms in a Mol which have it'''
    getter_type = RDPROP_GETTERS[prop_type]
    for atom in rdmol.GetAtoms():
        prop_val = getattr(atom, getter_type)(prop) # invoke type-appropriate getter on atom, with the name of the desired property
        
        if hasattr(prop_val, '__round__') and annotate_precision is not None: # only round on roundable objects, and only when
            prop_val = round(prop_val, annotate_precision)
        atom.SetProp('atomNote', str(prop_val)) # need to convert to string, as double is susceptible to float round display errors (shows all decimal places regardless of rounding)
    
@optional_in_place
def clear_atom_props(rdmol : Mol) -> None:
    '''Wipe properties of all atoms in a molecule'''
    for atom in rdmol.GetAtoms():
        for prop_name in atom.GetPropNames():
            atom.ClearProp(prop_name)