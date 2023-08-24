'''Utilities for finding, creating, and manipulating bonding sites (AKA "ports")'''

from typing import Generator
from rdkit import Chem

from . import bonding
from ..rdtypes import RDAtom, RDMol
from ..labeling import atomwise, molwise

from ...genutils.decorators.functional import optional_in_place
from ...genutils.iteration import iter_len


# CUSTOM EXCEPTIONS
class MolPortError(Exception):
    '''Raised when port-related errors as encountered'''
    pass

# DEFINING PORTS
PORT = Chem.AtomFromSmarts('[*]') # cast as lambda to generate new atom on each call (prevents cross-mutability issues)

def is_port(atom : RDAtom) -> bool: # TODO : add check for one single-bond connection
    '''Shorthand for deciding if an atom is a bonding site ("port")''' 
    return atom.GetAtomicNum() == 0 # NOTE : can't check equality with PORT atom query (would match every atom)

# QUERYING PORTS
def get_ports(rdmol : RDMol) -> Generator[RDAtom, None, None]:
    '''Generate all atoms which are ports'''
    for atom in rdmol.GetAtoms():
        if is_port(atom):  # NOTE : can't check equality with PORT atom query (would match every atom)
            yield atom

def get_port_ids(rdmol : RDMol) -> list[int]:
    '''Get indices of all atoms which are ports'''
    return [atom.GetIdx() for atom in get_ports(rdmol)]

def get_num_ports(rdmol : RDMol) -> int:
    '''Counts the number of port atoms present in a Mol'''
    return iter_len(get_ports(rdmol)) # NOTE : purposely avoid use of "len" to permit general iterators

def get_single_port(rdmol : RDMol) -> RDAtom:
    '''Get the singular port of a Mol which contains only 1 port'''
    num_ports = get_num_ports(rdmol)
    if num_ports == 0:
        raise MolPortError("Can't return port for molecule which has no ports")
    
    if num_ports > 1:
        raise MolPortError("Molecule has multiple ports; choice of single port is ambiguous")
    
    return next(get_ports(rdmol)) # return if port count checks pass

neighbor_ports     = atomwise._neighbor_factory_by_condition(    condition=is_port)
has_neighbor_ports = atomwise._has_neighbor_factory_by_condition(condition=is_port)

# PORT MODIFICATION
@optional_in_place
def hydrogenate_rdmol_ports(rdmol : RDMol) -> None:
    '''Replace all port atoms with hydrogens'''
    for port_id in get_port_ids(rdmol):
        rdmol.GetAtomWithIdx(port_id).SetAtomicNum(1)

def splice_port(rdmol : RDMol, sat_group : RDMol) -> RDMol:
    '''Splice a capping group onto a monomer port with matching port designation'''
    # PHASE 1 : labelling and locating bonding ports
    rdmol, sat_group = molwise.assign_contiguous_atom_map_nums(rdmol, sat_group, in_place=False) # VITAL that this is done first to ensure map nums are unique, while also unchanged during combination
    sat_port = get_single_port(sat_group)
    port_desig = sat_port.GetIsotope()
    sat_port_map_num = sat_port.GetAtomMapNum()

    for mono_port in get_ports(rdmol):
        if mono_port.GetIsotope() == port_desig:
            mono_port_map_num = mono_port.GetAtomMapNum()
            break # find first available port
    else:
        raise MolPortError(f'No compatible ports of designation "{port_desig}" found in molecule to be saturated')

    # PHASE 2a : preparing molecules for reaction
    combomol = Chem.CombineMols(rdmol, sat_group) # combine mols into single manipulable object
    combomol = Chem.RWMol(combomol) # make combo mutable for reaction manipulation

    # PHASE 2b : locating port-neighbor atoms to bond together
    bond_atom_ids = []
    for port_id in molwise.atom_ids_by_map_nums(combomol, sat_port_map_num, mono_port_map_num):
        port = combomol.GetAtomWithIdx(port_id)
        atom_to_bond = port.GetNeighbors()[0]
        bond_atom_ids.append(atom_to_bond.GetIdx())

    product = bonding.increase_bond_order(combomol, *bond_atom_ids, prioritize_unlabelled_ports=False) # must be set to False to ensure correct bonding location
    return RDMol(product) # convert back to immutable Mol and return

def saturate_ports(rdmol : RDMol, sat_group : RDMol) -> RDMol:
    '''Fill all possible ports which match a capping group with a particular port designation'''
    product = rdmol # initialize product as starting molecule (if no compatible ports are found, is returned unchanged)

    targ_desig = get_single_port(sat_group).GetIsotope()
    for port in get_ports(product):
        if port.GetIsotope() == targ_desig:
            product = splice_port(product, sat_group) # iteratively splice onto every compatible port

    return product 