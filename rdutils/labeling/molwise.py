'''For reading, writing, and clearing labels from all Atoms and/or Bonds in an RDKit molecule'''

from typing import Generator, Iterable, Optional, Union

from ..rdtypes import RDMol
from ...genutils.decorators.functional import optional_in_place


# READING FUNCTIONS
def get_isotopes(rdmol : RDMol, unique : bool=True) -> Union[set[int], list[int]]:
    '''Return all isotope IDs present in an RDMol. Can optionally return only the unique IDs'''
    isotope_ids = [atom.GetIsotope() for atom in rdmol.GetAtoms()]

    if unique:
        return set(isotope_ids)
    return isotope_ids

def ordered_map_nums(rdmol : RDMol) -> list[int]:
    '''Get assigned atom map numbers, in the same order as the internal RDMol atom IDs'''
    return [atom.GetAtomMapNum() for atom in rdmol.GetAtoms()]

def map_nums_by_atom_ids(rdmol : RDMol, *query_atom_ids : list[int]) -> Generator[Optional[int], None, None]: # TODO : generalize this to handle case where multiple atoms have the same map num
    '''Get assigned atom map numbers for a collection of atom ids, in the same order as the internal RDMol atom IDs'''
    for atom_id in query_atom_ids:
        yield(rdmol.GetAtomWithIdx(atom_id).GetAtomMapNum())

def atom_ids_by_map_nums(rdmol : RDMol, *query_map_nums : list[int]) -> Generator[Optional[int], None, None]: # TODO : generalize this to handle case where multiple atoms have the same map num
    '''Returns the first occurences of the atom IDs of any number of atoms, indexed by atom map number'''
    present_map_nums : list[int] = ordered_map_nums(rdmol)
    for map_num in query_map_nums:
        try:
            yield present_map_nums.index(map_num)
        except ValueError: # if the provided map number is not found, return NoneType
            yield None

# CHECKING FUNCTIONS
def has_fully_mapped_atoms(rdmol : RDMol) -> bool:
    '''Check whether an RDMol has a map number explicitly assigned to each member RDAtom'''
    for atom in rdmol.GetAtoms():
        if atom.GetAtomMapNum() == 0:
            return False
    else:
        return True
    
def has_uniquely_mapped_atoms(rdmol : RDMol) -> bool:
    '''Check whether an RDMol has distinct atom map numbers for each member RDAtom'''
    map_nums = ordered_map_nums(rdmol)
    return (len(map_nums) == len(set(map_nums))) # NOTE : not using rdmol.GetNumAtoms() as reference to avoid ambiguity with SMILES without explicit Hs


# WRITING FUNCTIONS
@optional_in_place    
def assign_ordered_atom_map_nums(rdmol : RDMol, start_from : int=1) -> None:
    '''Assigns atom's index as its atom map number for all atoms in an RDmol
    Can optionally specify what value to begin counting from (by default 1)'''
    for atom in rdmol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + start_from) # NOTE that starting from anything below 1 will cause an atom somewhere to be mapped to 0 (i.e. not mapped)

@optional_in_place    
def assign_atom_map_nums_from_ref(rdmol : RDMol, ref : Iterable[int]) -> None:
    '''Assigns atom map numbers using an external collection of values'''
    for atom, map_num in zip(rdmol.GetAtoms(), ref): # TODO : add some way to check that lengths match (may be generator-like)
        atom.SetAtomMapNum(map_num) 

# NOTE : this deliberately does NOT have an optional_in_place decorator (is implemented internally due to Iterable input)
def assign_contiguous_atom_map_nums(*rdmols : Iterable[RDMol], start_from : int=1, in_place : bool=False) -> Optional[list[RDMol]]: 
    '''Assign sequential numbering to a collection of molecules such that their map numbers span a contiguous range of integers
    Can optionally specify what value to begin counting from (by default 1)'''
    new_mols = []

    map_num_counter = start_from # initialize index tracker with starting index
    for rdmol in rdmols:
        new_mol = assign_ordered_atom_map_nums(rdmol, start_from=map_num_counter, in_place=in_place) # note : will be assigned NoneType if in-place
        if not in_place:
            new_mols.append(new_mol)
        map_num_counter += rdmol.GetNumAtoms()

    if new_mols:
        return new_mols


# CLEARING FUNCTIONS
@optional_in_place
def clear_atom_map_nums(rdmol : RDMol) -> None:
    '''Removes atom map numbers from all atoms in an RDMol'''
    for atom in rdmol.GetAtoms():
        atom.SetAtomMapNum(0)

@optional_in_place
def clear_atom_isotopes(rdmol : RDMol) -> None:
    '''Removes isotope numbers from all atoms in an RDMol'''
    for atom in rdmol.GetAtoms():
        atom.SetIsotope(0)