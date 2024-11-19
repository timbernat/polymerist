'''Tools for replacing Ports and functional groups'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Optional

from rdkit import Chem
from rdkit.Chem import Mol, RWMol

from ._bonding import combined_rdmol
from .formation import increase_bond_order
from .portlib import get_linker_ids, get_single_port
from .identification import get_num_bondable_port_pairs, get_first_bondable_port_pair

from ..labeling import molwise
from ...genutils.decorators.functional import optional_in_place


@optional_in_place
def splice_atoms(rwmol : RWMol, atom_id_1 : Optional[int]=None, atom_id_2 : Optional[int]=None, flavor_pair : tuple[Optional[int], Optional[int]]=(None, None)) -> None:
    '''Completely combine a single pair of ports on two target atoms to create a bond of desired order'''
    port_1, port_2 = get_first_bondable_port_pair(rwmol, atom_id_1=atom_id_1, atom_id_2=atom_id_2, flavor_pair=flavor_pair) # TODO : reimplement this to be slightly less redundant (maybe add option to specify bond order ahead of time in increase_bond_order()?)
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

def saturate_ports(rdmol : Mol, cap : Mol=Chem.MolFromSmiles('[*]-[H]'), flavor_to_saturate : int=0) -> None:
    '''Takes an RDKit Mol and another "cap" molecule (by default just hydrogen) and caps all valid ports (with the specified flavor) on the target Mol with the cap group'''
    flavor_pair : tuple[int, int] = (flavor_to_saturate, get_single_port(cap).flavor) # will raise exception if cap has anything other than 1 port
    
    rwmol = combined_rdmol(rdmol, cap, editable=True) # create initial combination to test if bonding is possible
    num_bonds_formable = get_num_bondable_port_pairs(rwmol, flavor_pair=flavor_pair)
    if num_bonds_formable == 0:
        return rdmol # special case, if no bondable ports exist to begin with, return the original molecule
    
    # implicit else if initial return in null case is not hit
    while num_bonds_formable > 0: # not greater than 1, since initial combination is already done outside loop
        splice_atoms(rwmol, flavor_pair=flavor_pair, in_place=True)
        if num_bonds_formable > 1: # need -1 to exclude final addition
            rwmol = combined_rdmol(rwmol, cap, editable=True) # add another cap group on all but the final iteration. NOTE : must be added one-at-a-time to preclude caps from bonding to each other
        num_bonds_formable = get_num_bondable_port_pairs(rwmol, flavor_pair=flavor_pair) # update count of available bonds (this may change as new bonds are added)
    molwise.assign_ordered_atom_map_nums(rwmol, in_place=True) # ensure map numbers are ordered and minial
    
    return Chem.Mol(rwmol) # revert to "regular" Mol from RWMol

@optional_in_place # temporarily placed here for backwards-compatibility reasons
def hydrogenate_rdmol_ports(rdmol : Mol) -> None:
    '''Replace all port atoms with hydrogens'''
    # return Chem.ReplaceSubstructs(rdmol, Chem.MolFromSmarts('[#0]'), Chem.MolFromSmarts('[#1]'), replaceAll=True)[0] # changes order of atom ids, undesiredly
    for port_id in get_linker_ids(rdmol):
        rdmol.GetAtomWithIdx(port_id).SetAtomicNum(1)