'''For reading, writing, and clearing labels from RDKit Atoms, Bonds, and Mols'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Generator, Iterable, Optional, Union
from ast import literal_eval

from rdkit.Chem.rdchem import Mol, Bond
from rdkit.Chem.rdmolfiles import MolToSmiles

from ..genutils.decorators.functional import optional_in_place


# CHECKING FUNCTIONS
def has_fully_mapped_atoms(rdmol : Mol) -> bool:
    '''Check whether an RDKit Mol has a map number explicitly assigned to each member Atom'''
    for atom in rdmol.GetAtoms():
        if atom.GetAtomMapNum() == 0:
            return False
    else:
        return True
    
def has_uniquely_mapped_atoms(rdmol : Mol, skip_unmapped : bool=False) -> bool: 
    '''
    Check whether an RDKit Mol has distinct atom map numbers for each member Atom
    If skip_unmapped=False (default), will check map numbers on ALL atoms;
    If skip_unmapped=True, however, will only check uniqueness of NONZERO map number (i.e. explicitly-mapped) atoms
    '''
    seen_map_numbers = set()
    for atom in rdmol.GetAtoms():
        map_num = atom.GetAtomMapNum()
        if (map_num == 0) and skip_unmapped:
            continue
        
        if map_num in seen_map_numbers:
            return False
        seen_map_numbers.add(map_num)
    else:
        return True
    
# READING FUNCTIONS
def map_numbers_by_atom_idxs(rdmol : Mol, *atom_idxs : list[int]) -> Generator[Optional[int], None, None]: # TODO : generalize this to handle case where multiple atoms have the same map num
    '''Get assigned atom map numbers for a collection of atom ids, in the same order as the internal RDKit Mol atom IDs'''
    for atom_idx in atom_idxs:
        yield(rdmol.GetAtomWithIdx(atom_idx).GetAtomMapNum())

def atom_idxs_by_map_numbers(rdmol : Mol, *map_numbers : list[int]) -> Generator[Optional[int], None, None]: # TODO : generalize this to handle case where multiple atoms have the same map num
    '''Returns the first occurences of the atom IDs of any number of atoms, indexed by atom map number'''
    present_map_nums : list[int] = [atom.GetAtomMapNum() for atom in rdmol.GetAtoms()]
    for map_num in map_numbers:
        try:
            yield present_map_nums.index(map_num)
        except ValueError: # if the provided map number is not found, return NoneType
            yield None
            
def get_bond_by_map_num_pair(rdmol : Mol, map_num_pair : tuple[int, int], as_bond : bool=True) -> Optional[Union[int, Bond]]:
    '''
    Get the bond spanning a pair of atoms with given pair of atom map numbers
    Returns the RDkit.Bond object if as_bond=True, and the index of the bond if as_bond=False
    
    If no bond exists between the atoms, will return None regardless of the value of "as_bond"
    '''
    bond = rdmol.GetBondBetweenAtoms(*atom_idxs_by_map_numbers(rdmol, *map_num_pair))
    if (not as_bond) and (bond is not None):
        return bond.GetIdx()
    return bond # returns bond or, implicitly, NoneType if no bond is found

# WRITING FUNCTIONS
@optional_in_place    
def assign_atom_map_nums_by_ids(rdmol : Mol, map_nums_by_ids : dict[int, int]) -> None:
    '''
    Assigns atom map numbers to Atoms in a Mol as identified by a mapping of atom indices to map numbers
    
    Parameters
    ----------
    rdmol : Chem.Mol
        An RDKit molecule to assign labels to
    map_nums_by_ids : dict[int, int]
        A dict keyed by atom index and mapping to the corresponding desired atom map numbers
    '''
    for atom_idx, map_num in map_nums_by_ids.items():
        rdmol.GetAtomWithIdx(atom_idx).SetAtomMapNum(map_num) 

@optional_in_place    
def assign_ordered_atom_map_nums(rdmol : Mol, start_from : int=1) -> None:
    '''Assigns atom's index as its atom map number for all atoms in an RDmol
    Can optionally specify what value to begin counting from (by default 1)'''
    assign_atom_map_nums_by_ids(
        rdmol,
        map_nums_by_ids={
            atom.GetIdx() : map_num
                for map_num, atom in enumerate(rdmol.GetAtoms(), start=start_from)
        },
        in_place=True, # assign_atom_map_nums_by_ids() already makes optional copies, so there's no need to make a second-order copy
    ) 

@optional_in_place
def relabel_map_nums(rdmol : Mol, relabeling : dict[int, int]) -> None:
    '''Applies a relabelling of atom map numbers to a subset of mapped atoms (described by a dict of old_map_num : new_map_num)'''
    assign_atom_map_nums_by_ids(
        rdmol,
        map_nums_by_ids={ # recast keys from current atom map nums to current atom ids (if even present)
            atom_idx : new_map_num
                for atom_idx, new_map_num in zip(atom_idxs_by_map_numbers(rdmol, *relabeling.keys()), relabeling.values())
                    if atom_idx is not None # TOSELF : consider adding check for duplicate remapping?
        },
        in_place=True, # assign_atom_map_nums_by_ids() already makes optional copies, so there's no need to make a second-order copy
    ) 
    
def mol_to_smiles_and_atom_permutation(mol: Mol, *args, **kwargs) -> tuple[str, list[int]]:
    '''
    Convert RDKit Mol to SMILES string (with any SMILES writer parameters passed to as args/kwargs)
    AND return permutation list which maps atoms in the written SMILES to their order in the passed Mol

    Useful when preserving Mol atom order in SMILES is necessary (not true in general)
    
    Parameters
    ----------
    mol : Chem.Mol
        The RDKit Mol object
    *args, **kwargs
        Additional arguments passed to the SMILES writer
        
    Returns
    -------
    smiles : str
        The resulting SMILES string
    atom_perm_inv : list[int]
        List representation of the permutation that restores atom order
        
        E.g. the following call, using "smiles" returned above, will in general scramble atom order:
        >>> from rdkit.Chem.rdmolfiles import MolToSmiles
        >>> mol = MolFromSmiles(smiles) # atom order doesn't match that of exporting Mol
        
        However, atom order can easily be restored by calling:
        >>> from rdkit.Chem.rdmolops import RenumberAtoms
        >>> mol = RenumberAtoms(mol, atom_perm_inv)

    tuple[str, list[int]]
        A tuple containing the SMILES string and the atom permutation list.
    '''
    smiles = MolToSmiles(mol, *args, **kwargs)
    # see RDKit "Magic" prop docs for more detail on this: https://www.rdkit.org/docs/RDKit_Book.html#romol-mol-in-python
    atom_perm : list[int] = literal_eval(mol.GetProp('_smilesAtomOutputOrder'))
    atom_perm_inv = sorted(range(mol.GetNumAtoms()), key=lambda i : atom_perm[i])
    
    return smiles, atom_perm_inv

# NOTE : this deliberately does NOT have an optional_in_place decorator (is implemented internally due to Iterable input)
def assign_contiguous_atom_map_nums(*rdmols : Iterable[Mol], start_from : int=1, in_place : bool=False) -> Optional[list[Mol]]: 
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
def clear_atom_map_nums(rdmol : Mol) -> None:
    '''Removes atom map numbers from all atoms in an RDKit Mol'''
    for atom in rdmol.GetAtoms():
        atom.SetAtomMapNum(0)

@optional_in_place
def clear_atom_isotopes(rdmol : Mol) -> None:
    '''Removes isotope numbers from all atoms in an RDKit Mol'''
    for atom in rdmol.GetAtoms():
        atom.SetIsotope(0)