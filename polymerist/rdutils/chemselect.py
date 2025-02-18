'''Utilities for conditional selection of chemical objects, such as atoms and bonds, from RDKit molecules'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Callable, Concatenate, Container, Union
from operator import (
    xor,
    xor as logical_xor, # alias for consistency
    or_  as logical_or,
    and_ as logical_and,
)
from rdkit import Chem


# CHEMICAL OBJECT TYPEHINTS
AtomCondition = Callable[Concatenate[Chem.Atom, ...], bool]
BondCondition = Callable[Concatenate[Chem.Bond, ...], bool]

AtomCollections = Union[set[int], set[Chem.Atom]]
BondCollections = Union[set[int], set[Chem.Bond], set[tuple[int, int]], set[tuple[Chem.Atom, Chem.Atom]]]

# CONDITIONAL SELECTION FUNCTIONS
def atoms_by_condition(
        mol : Chem.Mol,
        condition : AtomCondition=lambda atom : True,
        as_indices : bool=False,
        negate : bool=False,
    ) -> AtomCollections:
    '''
    Select a subset of atoms in a Mol based on a condition
    
    Parameters
    ----------
    mol : Chem.Mol
        An RDKit molecule object
    condition : Callable[[Chem.Atom], bool], default lambda atom : True
        Condition on atoms which returns bool; 
        Always returns True if unset
    as_indices : bool, default False
        Whether to return results as their indices (default) or as Atom objects
    negate : bool, default False
        Whether to invert the condition provided (by default False)
    
    Returns
    -------
    selected_atoms : Union[set[int], set[Chem.Atom]]
        A set of the atoms meeting the chosen condition
    '''
    return set(
        atom.GetIdx() if as_indices else atom
            for atom in mol.GetAtoms()
                if xor(condition(atom), negate)
    )

def bonds_by_condition(
        mol : Chem.Mol,
        condition : BondCondition=lambda bond : True,
        as_indices : bool=True,
        as_pairs : bool=True,
        negate : bool=False,
    ) -> BondCollections:
    '''
    Select a subset of bonds in a Mol based on a condition
    
    Parameters
    ----------
    mol : Chem.Mol
        An RDKit molecule object
    condition : Callable[[Chem.Bond], bool], default lambda bond : True
        Condition on bonds which returns bool; 
        Always returns True if unset
    as_indices : bool, default True
        Whether to return results as Bond objects or their indices (default)
    as_pairs : bool, default True
        Whether to return bonds as the pair of bondss they connect (default) or the bond itself
        Note that if as_pairs=True and as_indices=False, will return as pairs of Bonds objects
    negate : bool, default False
        Whether to invert the condition provided (by default False)
    
    Returns
    -------
    selected_bonds : Union[set[int], set[Chem.Bond]]
        A set of the bonds meeting the chosen condition
        Depending on flags set, bond will be represented as:
        * Bond indices
        * Bond objects
        * 2-tuples of Atom objects
        * 2-tuples of atom indices
    '''
    selected_bonds = set()
    for bond in mol.GetBonds():
        if xor(condition(bond), negate): # XOR
            if as_pairs:
                selected_bonds.add((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) if as_indices else (bond.GetBeginAtom(), bond.GetEndAtom()))
            else:
                selected_bonds.add(bond.GetIdx() if as_indices else bond)
                    
    return selected_bonds

# PREDEFINED CONDITIONS
def bond_condition_by_atom_condition_factory(atom_condition : AtomCondition, binary_operator : Callable[[bool, bool], bool]=logical_or) -> BondCondition:
    '''
    Dynamically define a bond condition based on an atom condition applied to the pair of atom a bond connects
    
    Evaluation over bond determined by a specified atom condition and a binary logical comparison made between the pair of atom condition evaluations
    By default, this binary condition is OR (i.e. the bond will evaluate True if either of its atoms meets the atom condition)
    '''
    def bond_condition(bond : Chem.Bond) -> bool:
        return binary_operator(atom_condition(bond.GetBeginAtom()), atom_condition(bond.GetEndAtom()))
    return bond_condition

# QUERIES BY SPECIFIC CONDITIONS
def get_mapped_atoms(mol : Chem.Mol, as_indices : bool=False) -> AtomCollections:
    '''Return all atoms (either as Atom objects or as indices) which have been assigned a nonzero atom map number'''
    return atoms_by_condition(
        mol,
        condition=lambda atom : atom.GetAtomMapNum() == 0,
        as_indices=as_indices,
        negate=True,
    )

def get_bonded_pairs(
        mol : Chem.Mol,
        *atom_idxs : Container[int],
        as_indices : bool=True,
        as_pairs : bool=True,
    ) -> BondCollections:
    '''Returns all bonds in a Mol which connect a pair of atoms whose indices both lie within the given atom indices'''
    return bonds_by_condition(
        mol,
        condition=bond_condition_by_atom_condition_factory(
            atom_condition=lambda atom : atom.GetIdx() in atom_idxs,
            binary_operator=logical_and,
        ),
        as_indices=as_indices,
        as_pairs=as_pairs,
        negate=False,
    )
    ...