'''Utilities for finding, creating, and manipulating bonding sites (AKA "ports")'''

from typing import Generator, Union

from itertools import product as cartesian_product
from collections import defaultdict
from dataclasses import dataclass

from rdkit import Chem

from ..reactions import bonding
from ..rdtypes import RDAtom, RDBond, RDMol
from ..labeling import atomwise, molwise

from ...genutils.decorators.functional import optional_in_place
from ...genutils.iteration import iter_len


# CUSTOM EXCEPTIONS
class MolPortError(Exception):
    '''Raised when port-related errors as encountered'''
    pass


# DEFINING A PORT
PORT = Chem.AtomFromSmarts('[*]') # cast as lambda to generate new atom on each call (prevents cross-mutability issues)

def is_port(atom : RDAtom) -> bool: # TODO : add check for one single-bond connection
    '''Shorthand for deciding if an atom is a bonding site ("port")''' 
    return (
        atom.GetAtomicNum() == 0 # check non-element atom type; NOTE : can't check equality with PORT atom query (would match every atom)
        and atomwise.get_num_bonds(atom) == 1 # ensure atom has on a single outgoign bond; NOTe : can't check ExplicitValence, since this might vary with bond order
    )


# QUERYING PORT(S)
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


# QUERYING ATOMS FOR NEIGHBOR PORTS
get_neighbor_ports = atomwise._get_neighbor_factory_by_condition(condition=is_port)
has_neighbor_ports = atomwise._has_neighbor_factory_by_condition(condition=is_port)

# QUERYING INFO ABOUT PORT(S)
@dataclass
class PortInfo:
    '''Useful information about a port regarding bonding'''
    port : RDAtom
    desig : int        # designation of port, indicating bonding compatibility
    bh_atom : RDAtom   # bridgehead atom, only neighboring "real" atom
    inc_bond : RDBond  # bond type of the single incoming bond

def get_port_info(atom : RDAtom) -> PortInfo:
    '''Get useful info container about a port'''
    assert(is_port(atom)) # necessary to ensure bridgehead atom and incoming bond are well-defined (single-neighbor only)

    return PortInfo(
        port=atom,
        desig=atom.GetIsotope(),
        bh_atom=atom.GetNeighbors()[0], 
        inc_bond=atom.GetBonds()[0]
    )

def get_mol_ports_dict(rdmol : RDMol) -> dict[int, list[PortInfo]]:
    '''Returns dict of all port info for a Mol, keyed by port designation'''
    ports_dict = defaultdict(list)
    for port in get_ports(rdmol):
        port_info = get_port_info(port)
        ports_dict[port_info.desig].append(port_info)

    return ports_dict

def get_atom_ports_dict(atom : RDAtom) -> dict[int, list[PortInfo]]:
    '''Returns dict of port info for all ports neighboring an atom, keyed by port designation'''
    ports_dict = defaultdict(list)
    for port in get_neighbor_ports(atom):
        port_info = get_port_info(port)
        ports_dict[port_info.desig].append(port_info)

    return ports_dict


# IDENTIFYING PORT BONDING INFORMATION
def ports_are_bondable(port_1 : RDAtom, port_2 : RDAtom) -> bool:
    '''Determine if two port atoms can be combined into a bond'''
    port_info_1 = get_port_info(port_1)
    port_info_2 = get_port_info(port_2)

    return (
        port_info_1.desig == port_info_2.desig # ensure that port designations are compatible... TODO : consider adding wild-like behavior for desig 0
        and port_info_1.inc_bond.GetBondType() ==  port_info_2.inc_bond.GetBondType() # ...and that bond types match
    )

def enumerate_bondable_atom_pairs(rdmol_1 : RDMol, rdmol_2 : RDMol, asAtoms : bool=True
                                  ) -> Union[list[tuple[tuple[int, int], tuple[int, int]]], list[tuple[tuple[RDAtom, RDAtom], tuple[RDAtom, RDAtom]]]]:
    '''Get all pairs of atoms between two Mols which have compatible neighboring ports
    Returns a dict with keys containing the bondable atoms and values containing the corresponding bond ports'''
    ports_dict_1 = get_mol_ports_dict(rdmol_1)
    ports_dict_2 = get_mol_ports_dict(rdmol_2)

    pairs_list = [
        ((port_info_1.bh_atom, port_info_2.bh_atom), (port_info_1.port, port_info_2.port))
            for mutual_desig in (ports_dict_1.keys() | ports_dict_2.keys()) # enumerate over pairs of ports which match designation
                for port_info_1, port_info_2 in cartesian_product(ports_dict_1[mutual_desig], ports_dict_2[mutual_desig]) # iterate over all pairs with matching designation
                    if ports_are_bondable(port_info_1.port, port_info_2.port)
    ]

    if not asAtoms:
        return [
            (tuple(atom.GetIdx() for atom in atom_pair), tuple(port.GetIdx() for port in port_pair))
                for atom_pair, port_pair in pairs_list
        ]
    return pairs_list

def num_bonds_formable(rdmol_1 : RDMol, rdmol_2 : RDMol) -> int:
    '''Return max number of new bonds that can possibly be formed between a pair of Mols'''
    return len(enumerate_bondable_atom_pairs(rdmol_1, rdmol_2))

def are_bondable(rdmol_1 : RDMol, rdmol_2 : RDMol) -> bool:
    '''Determine if any new bonds can be formed between a pair of Mols'''
    return num_bonds_formable(rdmol_1, rdmol_2) > 0
    # return any(enumerate_bondable_atom_pairs(rdmol_1, rdmol_2))


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