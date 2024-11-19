'''Tools for creating new bonds from free Ports in RDKit Mols'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import logging
LOGGER = logging.getLogger(__name__)

from typing import Optional
from rdkit.Chem.rdchem import BondType, RWMol

from .dissolution import _decrease_bond_order
from .identification import get_first_bondable_port_pair

from ..rdkdraw import clear_highlights
from ...genutils.decorators.functional import optional_in_place
from ...smileslib.primitives import RDKIT_BONDS_BY_BONDTYPE


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
        rwmol.ReplaceBond(curr_bond.GetIdx(), RDKIT_BONDS_BY_BONDTYPE[new_bond_type], preserveProps=True)

@optional_in_place
def increase_bond_order(rwmol : RWMol, atom_id_1 : int, atom_id_2 : int, flavor_pair : tuple[Optional[int], Optional[int]]=(None, None)) -> None: # TODO : add specificity to flavor selection
    '''Exchange two ports for a bond of one higher order in a modifiable RWMol. Can optionally specify a port flavor for greater selectivity'''
    port_pair = get_first_bondable_port_pair(rwmol, atom_id_1=atom_id_1, atom_id_2=atom_id_2, flavor_pair=flavor_pair) # raises MolPortError if none exists
    _increase_bond_order(rwmol, atom_id_1, atom_id_2, in_place=True) # Up-convert bond between target atoms; CRITICAL that this be done FIRST to avoid index shifts if linker atoms need to be removed
    
    # down-convert bonds to ports, removing linkers if necessary
    for port in port_pair:
        is_removed = (port.bond.GetBondType() == BondType.SINGLE) # check bond order before down-conversion
        _decrease_bond_order(rwmol, port.linker.GetIdx(), port.bridgehead.GetIdx(), in_place=True)
        if is_removed: # if the port bond vanishes when down-converting bond order, the linker atom must also be deleted
            rwmol.RemoveAtom(port.linker.GetIdx())

    clear_highlights(rwmol) # !NOTE! : this is necessary to allow for correct display in Jupyter if atoms are removed
