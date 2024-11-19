'''Tools for breaking bonds in RDKit Mols and assigning new Ports'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import logging
LOGGER = logging.getLogger(__name__)

from typing import Optional
from rdkit import Chem
from rdkit.Chem.rdchem import BondType, RWMol

from ..rderrors import BondOrderModificationError
from ..labeling.bondwise import are_bonded_atoms
from ...genutils.decorators.functional import optional_in_place


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
        new_linker = Chem.AtomFromSmarts('[#0]') # TODO : see if this can be made into a SMILES mol instead?
        new_linker.SetIsotope(flavor)
        new_port_id = rwmol.AddAtom(new_linker)# insert new port into molecule, taking note of index (TOSELF : ensure that this inserts indices at END of existing ones, could cause unexpected modification if not)
        
        rwmol.AddBond(atom_id, new_port_id, order=BondType.SINGLE) # bond the atom to the new port
        # _increase_bond_order(rwmol, atom_id, new_port_id)

@optional_in_place
def dissolve_bond(rwmol : RWMol, atom_id_1 : int, atom_id_2 : int, new_flavor_pair : Optional[tuple[int, int]]=None) -> None:
    '''Completely decompose a bond between two atoms, filling in ports with the chosen flavor''' 
    while rwmol.GetBondBetweenAtoms(atom_id_1, atom_id_2) is not None:
        decrease_bond_order(rwmol, atom_id_1, atom_id_2, new_flavor_pair=new_flavor_pair, in_place=True)
