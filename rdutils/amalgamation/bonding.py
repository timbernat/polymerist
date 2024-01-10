'''Tools for making and breaking bonds, with correct conversions of ports'''

import logging
LOGGER = logging.getLogger(__name__)

from typing import Iterable, Optional, Union
from functools import reduce
from collections import Counter
from IPython.display import display # for Jupyter display support

from rdkit import Chem
from rdkit.Chem.rdchem import BondType

from . import portlib
from ..smileslib.primitives import BOND_SMILES_BY_ORDER

from ..rdtypes import RWMol, RDMol, RDAtom
from ..rdkdraw import clear_highlights
from ..labeling import molwise

from ...genutils.maths.combinatorics.sequences import int_complement
from ...genutils.decorators.functional import optional_in_place


# BOND REFERENCE
class BondOrderModificationError(Exception):
    '''Raised when an invalid RDBond modification is attempted'''
    pass

def are_bonded_atoms(rdmol : RDMol, atom_id_1 : int, atom_id_2 : int) -> bool:
    '''Check if pair of atoms in an RDMol have a bond between then'''
    return (rdmol.GetBondBetweenAtoms(atom_id_1, atom_id_2) is not None)

def combined_rdmol(*rdmols : Iterable[RDMol], assign_map_nums : bool=True, editable : bool=True) -> Union[RDMol, RWMol]:
    '''Merge two RDMols into a single molecule with contiguous, non-overlapping atom map numbers'''
    if assign_map_nums: # assign contiguous, unique atom map numbers
        rdmols = molwise.assign_contiguous_atom_map_nums(*rdmols, in_place=False) 
    
    combo = reduce(Chem.CombineMols, rdmols) # combine into single Mol object to allow for bonding

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
        rwmol.ReplaceBond(curr_bond.GetIdx(), BOND_SMILES_BY_ORDER[new_bond_type], preserveProps=True)


# SINGLE BOND-ORDER CHANGES 
@optional_in_place
def decrease_bond_order(rwmol : RWMol, atom_id_1 : int, atom_id_2 : int, new_flavor_pair : Optional[tuple[int, int]]=None) -> None: 
    '''Lower the order of a bond between two atoms and insert two new ports in its place, raising Exception if no bond exists'''
    atom_ids = (atom_id_1, atom_id_2)
    if new_flavor_pair is None:
        new_flavor_pair = (0, 0)
    assert(len(new_flavor_pair) == 2)
    
    _decrease_bond_order(rwmol, atom_id_1, atom_id_2, in_place=True) # NOTE : must explicitly be called in-place to ensure correct top-level behavior, since this function is also decorated
    # free_isotope_labels = int_complement(get_isotopes(rwmol, unique=True), bounded=False) # generate unused isotope labels
    # add new ports for broken bond
    for atom_id, flavor in zip(atom_ids, new_flavor_pair):
        new_linker = Chem.AtomFromSmarts('[#0]') 
        new_linker.SetIsotope(flavor)
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
def dissolve_bond(rwmol : RWMol, atom_id_1 : int, atom_id_2 : int, new_flavor_pair : int=0) -> None:
    '''Completely decompose a bond between two atoms, filling in ports with the chosen flavor''' 
    while rwmol.GetBondBetweenAtoms(atom_id_1, atom_id_2) is not None:
        decrease_bond_order(rwmol, atom_id_1, atom_id_2, new_flavor_pair=new_flavor_pair, in_place=True)

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
    # return Chem.ReplaceSubstructs(rdmol, Chem.MolFromSmarts('[#0]'), Chem.MolFromSmarts('[#1]'), replaceAll=True)[0]
    for port_id in portlib.get_linker_ids(rdmol):
        rdmol.GetAtomWithIdx(port_id).SetAtomicNum(1)


# BOND SWAPS
def _is_valid_bond_derangement(bond_derangement : dict[int, tuple[int, int]]) -> bool:
    '''Determine whether an interatomic bond remapping describes a valid derangement'''
    # 1) check that each swap maps to a new element (i.e. no identity swaps)
    for begin_map_num, (curr_end_map_num, targ_end_map_num) in bond_derangement.items():
        if curr_end_map_num == targ_end_map_num:
            LOGGER.warn(f'Swap defined for initial index {begin_map_num} maps back to current partner ({curr_end_map_num} -> {targ_end_map_num})')
            return False
        
    # 2) check bijection (i.e. terminal atom remappings form a closed multiset)
    curr_end_counts, targ_end_counts = [Counter(i) for i in zip(*bond_derangement.values())] #  multisets are permissible for when multiple current/target bonds connect to the same atom 
    if curr_end_counts != targ_end_counts:
        LOGGER.warn('Bond derangement does not define a 1-1 correspondence between elements in the multiset')
        return False

    return True # only return if all aabove checks pass

@optional_in_place
def swap_bonds(rwmol : RWMol, bond_derangement : dict[int, tuple[int, int]], show_steps : bool=False) -> Optional[RWMol]:
    '''
    Takes a modifiable Mol and a bond derangement dict and performs the requested bond swaps
    Derangement dict should have th following form:
        keys   : int             = corresponds to the beginning atom of a bond
        values : tuple[int, int] = corresponds to the current end atom map number and target end atom map number (in that order) 
    Modifiable Mol can contain multiple disconnected molecular components
    ''' 
    # TODO : check for complete atom map num assignment
    if not _is_valid_bond_derangement(bond_derangement):
        raise ValueError('Invalid interatomic bond derangement provided')

    # determine non-interfering port flavors for new bonds (preserves parity between permutation sets)
    available_port_flavors = int_complement(molwise.get_isotopes(rwmol), bounded=False) # ensures newly-created temporary ports don't clash with any existing ones
    flavor_pair = (next(available_port_flavors), next(available_port_flavors)) # grab first two available flavors
    portlib.Port.bondable_flavors.insert(flavor_pair) # temporarily register pair as bondable

    # break current bonds
    for begin_map_num, (curr_end_map_num, _) in bond_derangement.items():
        decrease_bond_order(
            rwmol,
            *molwise.atom_ids_by_map_nums(rwmol, begin_map_num, curr_end_map_num),
            new_flavor_pair=flavor_pair,
            in_place=True # must be done in-place to allow optional_in_place decoration
        )

        if show_steps:
            print(f'{begin_map_num} --x-> {curr_end_map_num}')
            display(rwmol)

    # form new bonds - must be done AFTER breakage to ensure all necessary ports exist
    for begin_map_num, (_, targ_end_map_num) in bond_derangement.items():
        increase_bond_order(
            rwmol,
            *molwise.atom_ids_by_map_nums(rwmol, begin_map_num, targ_end_map_num), 
            flavor_pair=flavor_pair,
            in_place=True # must be done in-place to allow optional_in_place decoration
        )

        if show_steps:
            print(f'{begin_map_num} ----> {targ_end_map_num}')
            display(rwmol)

    # deregister bondable pair
    portlib.Port.bondable_flavors.pop(flavor_pair) 