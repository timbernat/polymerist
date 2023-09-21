'''Tools for making and breaking bonds, with correct conversions of ports'''

from typing import Optional, Union

from rdkit import Chem
from rdkit.Chem.rdchem import BondType

from .portlib import get_bondable_port_pair_between_atoms, max_bondable_order_between_atoms
from ..rdtypes import RWMol, RDMol, RDAtom
from ..rdkdraw import clear_highlights
from ..labeling import molwise
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

def are_bonded_atoms(rdmol : RDMol, atom_id_1 : int, atom_id_2 : int) -> bool:
    '''Check if pair of atoms in an RDMol have a bond between then'''
    return (rdmol.GetBondBetweenAtoms(atom_id_1, atom_id_2) is not None)

def combine_rdmols(rdmol_1 : RDMol, rdmol_2 : RDMol, editable : bool=True) -> Union[RDMol, RWMol]:
    '''Merge two RDMols into a single molecule with contiguous, non-overlapping atom map numbers'''
    rdmol_1, rdmol_2 = molwise.assign_contiguous_atom_map_nums(rdmol_1, rdmol_2, in_place=False) 
    combo = Chem.CombineMols(rdmol_1, rdmol_2) # combine into single Mol object to allow for bonding

    if editable:
        return Chem.RWMol(combo) # make combined Mol modifiable
    return combo


# BASE BOND ORDER CHANGE FUNCTIONS
@optional_in_place
def _decrease_bond_order(rwmol : RWMol, atom_id_1 : int, atom_id_2 : int) -> RWMol: 
    '''Lower the order of a bond between two atoms, raising Expection if no bond exists
    DOES NOT ENSURE VALENCE OF BONDED ATOMS IS PRESERVED'''
    if not are_bonded_atoms(rwmol, atom_id_1, atom_id_2):
        raise BondOrderModificationError
    
    # determine expected bond type after order decrease (handle single-bond case, specifically) 
    curr_bond = rwmol.GetBondBetweenAtoms(atom_id_1, atom_id_2) # guaranteed to not be None by the bond_order_decreasable check at the start
    new_bond_type = BondType.values[curr_bond.GetBondTypeAsDouble() - 1] # with pre-existing bond, need to get the next order up by numeric lookup
    if new_bond_type == BondType.UNSPECIFIED:
        new_bond_type = None # explicitly set to NoneType if single bond is broken

    # remove existing bond; not single bond, replace with bond of new type
    rwmol.RemoveBond(atom_id_1, atom_id_2) # create new bond or specified order
    if new_bond_type is not None:
        rwmol.AddBond(atom_id_1, atom_id_2, order=new_bond_type)

@optional_in_place
def _increase_bond_order(rwmol : RWMol, atom_id_1 : int, atom_id_2 : int) -> None:
    '''Raise the order of a bond between two atoms, creating a new single-bond if none exists
    DOES NOT ENSURE VALENCE OF BONDED ATOMS IS PRESERVED'''
    # determine expected bond type after order increase (handle single-bond removal, specifically) 
    curr_bond = rwmol.GetBondBetweenAtoms(atom_id_1, atom_id_2)
    if curr_bond is None:
        new_bond_type = BondType.SINGLE # with no pre-existing bond, simply add a single bond
    else: 
        new_bond_type = BondType.values[curr_bond.GetBondTypeAsDouble() + 1] # with pre-existing bond, need to get the next order up by numeric lookup
        rwmol.RemoveBond(atom_id_1, atom_id_2) # also remove the existing bond for new bond creation

    rwmol.AddBond(atom_id_1, atom_id_2, order=new_bond_type) # create new bond or specified order

@optional_in_place
def _increase_bond_order_alt(rwmol : RWMol, atom_id_1 : int, atom_id_2 : int) -> None:
    '''ALTERNATE IMPLEMENTATION - Raise the order of a bond between two atoms, creating a new single-bond if none exists
    DOES NOT ENSURE VALENCE OF BONDED ATOMS IS PRESERVED'''
    # determine expected bond type after order increase (handle single-bond removal, specifically) 
    curr_bond = rwmol.GetBondBetweenAtoms(atom_id_1, atom_id_2)
    if curr_bond is None:
        rwmol.AddBond(atom_id_1, atom_id_2, order=BondType.SINGLE) # create new bond or specified order
    else: 
        new_bond_type = BondType.values[curr_bond.GetBondTypeAsDouble() + 1] # with pre-existing bond, need to get the next order up by numeric lookup
        rwmol.ReplaceBond(curr_bond.GetIdx(), BONDS_BY_ORDER[new_bond_type], preserveProps=True)


# SINGLE BOND-ORDER CHANGES 
@optional_in_place
def decrease_bond_order(rwmol : RWMol, atom_id_1 : int, atom_id_2 : int, new_port_desig : int=0) -> None: 
    '''Lower the order of a bond between two atoms and insert two new ports in its place, raising Expection if no bond exists'''
    _decrease_bond_order(rwmol, atom_id_1, atom_id_2, in_place=True) # NOTE : must explicitly be called in-place to ensure correct top-level behavior, since this function is also decorated
    
    # free_isotope_labels = int_complement(get_isotopes(rwmol, unique=True), bounded=False) # generate unused isotope labels
    # add new ports for broken bond
    for atom_id in (atom_id_1, atom_id_2):
        new_linker = Chem.AtomFromSmarts('[#0]') 
        new_linker.SetIsotope(new_port_desig)
        new_port_id = rwmol.AddAtom(new_linker)# insert new port into molecule, taking note of index (TOSELF : ensure that this inserts indices at END of existing ones, could cause unexpected modification if not)
        rwmol.AddBond(atom_id, new_port_id, order=BondType.SINGLE) # bond the atom to the new port
        # _increase_bond_order(rwmol, atom_id, new_port_id)

@optional_in_place
def increase_bond_order(rwmol : RWMol, atom_id_1 : int, atom_id_2 : int, target_desig : Optional[int]=None) -> None: # TODO : add specificity to designation selection
    '''Exchange two ports for a bond of one higher order in a modifiable RWMol. Can optionally specify a port designation for greater selectivity'''
    port_pair = get_bondable_port_pair_between_atoms(rwmol, atom_id_1, atom_id_2, target_desig=target_desig) # raises MolPortError if none exists

    # Up-convert bond between target atoms
    _increase_bond_order(rwmol, atom_id_1, atom_id_2, in_place=True) # important that new bond formation be done FIRST, to avoid index shifts if linker atoms need to be removed
    
    # down-convert bonds to ports, removing linkers if necessary
    for port in port_pair:
        is_removed = (port.bond.GetBondType() == BondType.SINGLE) # check bond order before down-conversion
        _decrease_bond_order(rwmol, port.linker.GetIdx(), port.bridgehead.GetIdx(), in_place=True)
        if is_removed: # if the port bond vanishes when down-converting bond order, the linker atom must also be deleted
            rwmol.RemoveAtom(port.linker.GetIdx())

    clear_highlights(rwmol) # !NOTE! : this is necessary to allow for correct display in Jupyter if atoms are removed


# TOTAL BOND CHANGES
@optional_in_place
def dissolve_bond(rwmol : RWMol, atom_id_1 : int, atom_id_2 : int, new_port_desig : int=0) -> None:
    '''Completely decompose a bond between two atoms, filling in ports with the chosen designation''' 
    while rwmol.GetBondBetweenAtoms(atom_id_1, atom_id_2) is not None:
        decrease_bond_order(rwmol, atom_id_1, atom_id_2, new_port_desig=new_port_desig, in_place=True)

@optional_in_place
def splice_atoms(rwmol : RWMol, atom_id_1 : int, atom_id_2 : int, target_desig : Optional[int]=None) -> None:
    '''Completely combine all bondable ports between a pair of atoms into one bond'''
    for i in range(max_bondable_order_between_atoms(rwmol, atom_id_1, atom_id_2, target_desig=target_desig)):
        if i == 0: # record map numbers as invariant between bonding events (atom id won't necessarily be)
            atom_map_nums = [j for j in molwise.map_nums_by_atom_ids(rwmol, atom_id_1, atom_id_2)] # unpack as list to avoid values being "used up" upon iteration
        else:
            atom_id_1, atom_id_2 = molwise.atom_ids_by_map_nums(rwmol, *atom_map_nums) # reassign bonded atom IDs from invariant map nums cached on first bond formation

        increase_bond_order(rwmol, atom_id_1, atom_id_2, target_desig=target_desig, in_place=True)