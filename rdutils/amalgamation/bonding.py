'''Tools for making and breaking bonds, with correct conversions of ports'''

from typing import Optional
from rdkit import Chem, BondType

from . import portlib, _bonding
from ..rdtypes import RWMol
from ..rdkdraw import clear_highlights

from ...genutils.decorators.functional import optional_in_place


# SINGLE BOND-ORDER CHANGES 
@optional_in_place
def decrease_bond_order(rwmol : RWMol, *bond_atom_ids : list[int, int], new_port_desig : int=0) -> None: 
    '''Lower the order of a bond between two atoms and insert two new ports in its place, raising Expection if no bond exists'''
    _bonding._decrease_bond_order(rwmol, *bond_atom_ids, in_place=True) # NOTE : must explicitly be called in-place to ensure correct top-level behavior, since this function is also decorated
    
    # free_isotope_labels = int_complement(get_isotopes(rwmol, unique=True), bounded=False) # generate unused isotope labels
    # add new ports for broken bond
    for atom_id in bond_atom_ids:
        new_port = Chem.AtomFromSmarts('[#0]') 
        new_port.SetIsotope(new_port_desig)
        new_port_id = rwmol.AddAtom(new_port)# insert new port into molecule, taking note of index (TOSELF : ensure that this inserts indices at END of existing ones, could cause unexpected modification if not)
        rwmol.AddBond(atom_id, new_port_id) # bond the atom to the new port
        # _bonding._increase_bond_order(rwmol, atom_id, new_port_id)

@optional_in_place
def increase_bond_order(rwmol : RWMol, *bond_atom_ids : list[int, int], port_desig : Optional[int]=None) -> None: # TODO : add specificity to designation selection
    '''Exchange two ports for a bond of one higher order in a modifiable RWMol'''
    id1, id2 = bond_atom_ids
    port_pairs = portlib.get_bondable_port_pairs(
        portlib.get_ports_on_atom_at_idx(rwmol, id1),
        portlib.get_ports_on_atom_at_idx(rwmol, id2)
    )

    # obtain first valid port pair
    try:
        port_pair = next(port_pairs) # no specificity (takes first bondable pair)
        port_pair = sorted(port_pair, key=lambda port : port.linker.GetIdx(), reverse=True) # sort in-place to ensure highest-index linker is first (MATTER FOR ORDER OF ATOM REMOVAL)
    except StopIteration:
        raise portlib.MolPortError(f'No bondable ports exist between atoms {id1} and {id2}')

    # Up-convert bond between target atoms
    _bonding._increase_bond_order(rwmol, *bond_atom_ids, in_place=True) # important that new bond formation be done FIRST, to avoid index shifts if linker atoms need to be removed
    
    # down-convert bonds to ports, removing linkers if necessary
    for port in port_pair:
        is_removed = (port.bond.GetBondType() == BondType.SINGLE) # check bond order before down-conversion
        _bonding._decrease_bond_order(rwmol, port.linker.GetIdx(), port.bridgehead.GetIdx(), in_place=True)
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
def fuse_atoms(rwmol : RWMol, *bond_atom_ids : list[int, int], port_desig : Optional[int]=None) -> None:
    '''Completely combine all bondable ports between a pair of atoms into one bond'''
    raise NotImplemented