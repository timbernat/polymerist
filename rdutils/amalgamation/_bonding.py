'''Bare-boned utilities for creating and destroying bonds between RDMols'''

from rdkit.Chem import RWMol, BondType

from ..rdtypes import RDMol, RDAtom
from ...genutils.decorators.functional import optional_in_place


class BondOrderModificationError(Exception):
    '''Raised when an invalid RDBond modification is attempted'''
    pass

# DOWN-CONVERTING BONDS
def are_bonded_atoms(rdmol : RDMol, *atom_pair_ids : list[int, int]) -> bool:
    '''Check if pair of atoms in an RDMol have a bond between then'''
    return (rdmol.GetBondBetweenAtoms(*atom_pair_ids) is not None)

@optional_in_place
def _decrease_bond_order(rwmol : RWMol, *bond_atom_ids : list[int, int]) -> RWMol: 
    '''Lower the order of a bond between two atoms, raising Expection if no bond exists
    DOES NOT ENSURE VALENCE OF BONDED ATOMS IS PRESERVED'''
    if not are_bonded_atoms(rwmol, *bond_atom_ids):
        raise BondOrderModificationError
    
    # determine expected bond type after order decrease (handle single-bond case, specifically) 
    curr_bond = rwmol.GetBondBetweenAtoms(*bond_atom_ids) # guaranteed to not be None by the bond_order_decreasable check at the start
    new_bond_type = BondType.values[curr_bond.GetBondTypeAsDouble() - 1] # with pre-existing bond, need to get the next order up by numeric lookup
    if new_bond_type == BondType.UNSPECIFIED:
        new_bond_type = None # explicitly set to NoneType if single bond is broken

    # remove existing bond; not single bond, replace with bond of new type
    rwmol.RemoveBond(*bond_atom_ids) # create new bond or specified order
    if new_bond_type is not None:
        rwmol.AddBond(*bond_atom_ids, order=new_bond_type)


# UP-CONVERTING ATOMS
@optional_in_place
def _increase_bond_order(rwmol : RWMol, *bond_atom_ids : list[int, int]) -> None:
    '''Raise the order of a bond between two atoms, creating a new single-bond if none exists
    DOES NOT ENSURE VALENCE OF BONDED ATOMS IS PRESERVED'''
    # determine expected bond type after order increase (handle single-bond removal, specifically) 
    curr_bond = rwmol.GetBondBetweenAtoms(*bond_atom_ids)
    if curr_bond is None:
        new_bond_type = BondType.SINGLE # with no pre-existing bond, simply add a single bond
    else: 
        new_bond_type = BondType.values[curr_bond.GetBondTypeAsDouble() + 1] # with pre-existing bond, need to get the next order up by numeric lookup
        rwmol.RemoveBond(*bond_atom_ids) # also remove the existing bond for new bond creation

    rwmol.AddBond(*bond_atom_ids, order=new_bond_type) # create new bond or specified order
