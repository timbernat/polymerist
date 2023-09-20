'''Tools for making and breaking bonds, with correct conversions of ports'''

from typing import Optional

from rdkit import Chem
from rdkit.Chem.rdchem import BondType

from .portlib import get_bondable_port_pair_between_atoms
from ..rdtypes import RWMol, RDMol, RDAtom
from ..rdkdraw import clear_highlights
from ...genutils.decorators.functional import optional_in_place


# BOND REFERENCE
BONDS_BY_ORDER = { # conrete bond objects (since these can't be directly instantiated from bond order in Python)
    BondType.SINGLE    : Chem.BondFromSmiles('-'),
    BondType.DOUBLE    : Chem.BondFromSmiles('='),
    BondType.TRIPLE    : Chem.BondFromSmiles('#'),
    BondType.QUADRUPLE : Chem.BondFromSmiles('$'),
    BondType.AROMATIC  : Chem.BondFromSmiles(':'),
}

class BondOrderModificationError(Exception):
    '''Raised when an invalid RDBond modification is attempted'''
    pass

def are_bonded_atoms(rdmol : RDMol, *atom_pair_ids : list[int, int]) -> bool:
    '''Check if pair of atoms in an RDMol have a bond between then'''
    return (rdmol.GetBondBetweenAtoms(*atom_pair_ids) is not None)


# BASE BOND ORDER CHANGE FUNCTIONS
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

@optional_in_place
def _increase_bond_order_alt(rwmol : RWMol, *bond_atom_ids : list[int, int]) -> None:
    '''ALTERNATE IMPLEMENTATION - Raise the order of a bond between two atoms, creating a new single-bond if none exists
    DOES NOT ENSURE VALENCE OF BONDED ATOMS IS PRESERVED'''
    # determine expected bond type after order increase (handle single-bond removal, specifically) 
    curr_bond = rwmol.GetBondBetweenAtoms(*bond_atom_ids)
    if curr_bond is None:
        rwmol.AddBond(*bond_atom_ids, order=BondType.SINGLE) # create new bond or specified order
    else: 
        new_bond_type = BondType.values[curr_bond.GetBondTypeAsDouble() + 1] # with pre-existing bond, need to get the next order up by numeric lookup
        rwmol.ReplaceBond(curr_bond.GetIdx(), BONDS_BY_ORDER[new_bond_type], preserveProps=True)


# SINGLE BOND-ORDER CHANGES 
@optional_in_place
def decrease_bond_order(rwmol : RWMol, *bond_atom_ids : list[int, int], new_port_desig : int=0) -> None: 
    '''Lower the order of a bond between two atoms and insert two new ports in its place, raising Expection if no bond exists'''
    _decrease_bond_order(rwmol, *bond_atom_ids, in_place=True) # NOTE : must explicitly be called in-place to ensure correct top-level behavior, since this function is also decorated
    
    # free_isotope_labels = int_complement(get_isotopes(rwmol, unique=True), bounded=False) # generate unused isotope labels
    # add new ports for broken bond
    for atom_id in bond_atom_ids:
        new_linker = Chem.AtomFromSmarts('[#0]') 
        new_linker.SetIsotope(new_port_desig)
        new_port_id = rwmol.AddAtom(new_linker)# insert new port into molecule, taking note of index (TOSELF : ensure that this inserts indices at END of existing ones, could cause unexpected modification if not)
        rwmol.AddBond(atom_id, new_port_id, order=BondType.SINGLE) # bond the atom to the new port
        # _increase_bond_order(rwmol, atom_id, new_port_id)

@optional_in_place
def increase_bond_order(rwmol : RWMol, *bond_atom_ids : list[int, int], targ_port_desig : Optional[int]=None) -> None: # TODO : add specificity to designation selection
    '''Exchange two ports for a bond of one higher order in a modifiable RWMol. Can optionally specify a port designation for greater selectivity'''
    id1, id2 = bond_atom_ids
    port_pair = get_bondable_port_pair_between_atoms(rwmol, id1, id2, targ_port_desig=targ_port_desig) # raises MolPortError if none exists

    # Up-convert bond between target atoms
    _increase_bond_order(rwmol, *bond_atom_ids, in_place=True) # important that new bond formation be done FIRST, to avoid index shifts if linker atoms need to be removed
    
    # down-convert bonds to ports, removing linkers if necessary
    for port in port_pair:
        is_removed = (port.bond.GetBondType() == BondType.SINGLE) # check bond order before down-conversion
        _decrease_bond_order(rwmol, port.linker.GetIdx(), port.bridgehead.GetIdx(), in_place=True)
        if is_removed: # if the port bond vanishes when down-converting bond order, the linker atom must also be deleted
            rwmol.RemoveAtom(port.linker.GetIdx())

    clear_highlights(rwmol) # !NOTE! : this is necessary to allow for correct display in Jupyter if atoms are removed


# TOTAL BOND CHANGES
@optional_in_place
def dissolve_bond(rwmol : RWMol, *bond_atom_ids : list[int, int], new_port_desig : int=0) -> None:
    '''Completely decompose a bond between two atoms, filling in ports with the chosen designation''' 
    while rwmol.GetBondBetweenAtoms(*bond_atom_ids) is not None:
        decrease_bond_order(rwmol, *bond_atom_ids, new_port_desig=new_port_desig, in_place=True)

@optional_in_place
def splice_atoms(rwmol : RWMol, *bond_atom_ids : list[int, int], port_desig : Optional[int]=None) -> None:
    '''Completely combine all bondable ports between a pair of atoms into one bond'''
    raise NotImplemented