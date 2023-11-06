'''Tools for making and breaking bonds, with correct conversions of ports'''

from typing import Optional, Union
from IPython.display import display
from rdkit import Chem
from rdkit.Chem.rdchem import BondType

from . import portlib
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

def combined_rdmol(rdmol_1 : RDMol, rdmol_2 : RDMol, editable : bool=True) -> Union[RDMol, RWMol]:
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
def decrease_bond_order(rwmol : RWMol, atom_id_1 : int, atom_id_2 : int, new_port_flavor : int=0) -> None: 
    '''Lower the order of a bond between two atoms and insert two new ports in its place, raising Expection if no bond exists'''
    _decrease_bond_order(rwmol, atom_id_1, atom_id_2, in_place=True) # NOTE : must explicitly be called in-place to ensure correct top-level behavior, since this function is also decorated
    
    # free_isotope_labels = int_complement(get_isotopes(rwmol, unique=True), bounded=False) # generate unused isotope labels
    # add new ports for broken bond
    for atom_id in (atom_id_1, atom_id_2):
        new_linker = Chem.AtomFromSmarts('[#0]') 
        new_linker.SetIsotope(new_port_flavor)
        new_port_id = rwmol.AddAtom(new_linker)# insert new port into molecule, taking note of index (TOSELF : ensure that this inserts indices at END of existing ones, could cause unexpected modification if not)
        rwmol.AddBond(atom_id, new_port_id, order=BondType.SINGLE) # bond the atom to the new port
        # _increase_bond_order(rwmol, atom_id, new_port_id)

@optional_in_place
def increase_bond_order(rwmol : RWMol, atom_id_1 : int, atom_id_2 : int, flavor_pair : tuple[Optional[int], Optional[int]]=(None, None)) -> None: # TODO : add specificity to flavor selection
    '''Exchange two ports for a bond of one higher order in a modifiable RWMol. Can optionally specify a port flavor for greater selectivity'''
    port_pair = portlib.get_first_bondable_port_pair(rwmol, atom_id_1=atom_id_1, atom_id_2=atom_id_2, flavor_pair=flavor_pair) # raises MolPortError if none exists
    _increase_bond_order(rwmol, atom_id_1, atom_id_2, in_place=True) # Up-convert bond between target atoms; CRITICAL that this be done FIRST to avoid index shifts if linker atoms need to be removed
    
    # down-convert bonds to ports, removing linkers if necessary
    for port in port_pair:
        is_removed = (port.bond.GetBondType() == BondType.SINGLE) # check bond order before down-conversion
        _decrease_bond_order(rwmol, port.linker.GetIdx(), port.bridgehead.GetIdx(), in_place=True)
        if is_removed: # if the port bond vanishes when down-converting bond order, the linker atom must also be deleted
            rwmol.RemoveAtom(port.linker.GetIdx())

    clear_highlights(rwmol) # !NOTE! : this is necessary to allow for correct display in Jupyter if atoms are removed


# TOTAL BOND CHANGES
@optional_in_place
def dissolve_bond(rwmol : RWMol, atom_id_1 : int, atom_id_2 : int, new_port_flavor : int=0) -> None:
    '''Completely decompose a bond between two atoms, filling in ports with the chosen flavor''' 
    while rwmol.GetBondBetweenAtoms(atom_id_1, atom_id_2) is not None:
        decrease_bond_order(rwmol, atom_id_1, atom_id_2, new_port_flavor=new_port_flavor, in_place=True)

@optional_in_place
def splice_atoms(rwmol : RWMol, atom_id_1 : Optional[int]=None, atom_id_2 : Optional[int]=None, flavor_pair : tuple[Optional[int], Optional[int]]=(None, None)) -> None:
    '''Completely combine a single pair of ports on two target atoms to create a bond of desired order'''
    port_1, port_2 = portlib.get_first_bondable_port_pair(rwmol, atom_id_1=atom_id_1, atom_id_2=atom_id_2, flavor_pair=flavor_pair) # TODO : reimplement this to be slightly less redundant (maybe add option to specify bond order ahead of time in increase_bond_order()?)
    if atom_id_1 is None: # fill in atom indices to give greater specificity for single bond order upconversion
        atom_id_1 = port_1.bridgehead.GetIdx()
    if atom_id_2 is None: # fill in atom indices to give greater specificity for single bond order upconversion
        atom_id_2 = port_2.bridgehead.GetIdx()

    assert(port_1.bond.GetBondType() == port_2.bond.GetBondType()) # slightly redundant given pre-bonding checks, but is a helpful failsafe
    for i in range(int(port_1.bond.GetBondTypeAsDouble())): # bond the target atoms up to the degree of the desired port pair; types
        if i == 0: # record map numbers, as these are invariant between bonding events (unlike atom IDs)
            atom_map_nums = [j for j in molwise.map_nums_by_atom_ids(rwmol, atom_id_1, atom_id_2)] # unpack as list to avoid values being "used up" upon iteration
        else:
            atom_id_1, atom_id_2 = molwise.atom_ids_by_map_nums(rwmol, *atom_map_nums) # reassign bonded atom IDs from invariant map nums cached on first bond formation

        increase_bond_order(rwmol, atom_id_1, atom_id_2, flavor_pair=flavor_pair, in_place=True) 

def saturate_ports(rdmol : RDMol, cap : RDMol=Chem.MolFromSmarts('[#0]-[#1]'), flavor_to_saturate : int=0) -> None:
    '''Takes an RDMol and another "cap" molecule (by default just hydrogen) and caps all valid ports (with the specified flavor) on the target Mol with the cap group'''
    flavor_pair : tuple[int, int] = (flavor_to_saturate, portlib.get_single_port(cap).flavor) # will raise exception if cap has anything other than 1 port
    
    rwmol = combined_rdmol(rdmol, cap, editable=True) # create initial combination to test if bonding is possible
    num_bonds_formable = portlib.get_num_bondable_port_pairs(rwmol, flavor_pair=flavor_pair)
    if num_bonds_formable == 0:
        return rdmol # special case, if no bondable ports exist to begin with, return the original molecule
    
    # implicit else if initial return in null case is not hit
    while num_bonds_formable > 0: # not greater than 1, since initial combination is already done outside loop
        splice_atoms(rwmol, flavor_pair=flavor_pair, in_place=True)
        if num_bonds_formable > 1: # need -1 to exclude final addition
            rwmol = combined_rdmol(rwmol, cap, editable=True) # add another cap group on all but the final iteration. NOTE : must be added one-at-a-time to preclude caps from bonding to each other
        num_bonds_formable = portlib.get_num_bondable_port_pairs(rwmol, flavor_pair=flavor_pair) # update count of available bonds (this may change as new bonds are added)
    molwise.assign_ordered_atom_map_nums(rwmol, in_place=True) # ensure map numbers are ordered and minial
    
    return Chem.rdchem.Mol(rwmol) # revert to "regular" Mol from RWMol

@optional_in_place # temporarily placed here for backwards-compatibility reasons
def hydrogenate_rdmol_ports(rdmol : RDMol) -> None:
    '''Replace all port atoms with hydrogens'''
    for port_id in portlib.get_linker_ids(rdmol):
        rdmol.GetAtomWithIdx(port_id).SetAtomicNum(1)