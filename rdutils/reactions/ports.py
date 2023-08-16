'''Utilities for finding, creating, and manipulating bonding sites (AKA "ports")'''

from ..rdtypes import RDAtom, RDMol
from ..labeling import atomwise

from ...genutils.decorators.functional import optional_in_place


# PORT INFO QUERY
def is_port(atom : RDAtom) -> bool:
    '''Shorthand for deciding if an atom is a bonding site ("port")''' 
    return atom.GetAtomicNum() == 0

def num_ports(rdmol : RDMol) -> int:
    '''Counts the number of port atoms present in a Mol'''
    return sum(1 for atom in rdmol.GetAtoms() if is_port(atom))

def get_port_ids(rdmol : RDMol) -> list[int]:
    '''Get atom indices of port (i.e. wild *-type or undefined) atoms'''
    return [
        atom.GetIdx()
            for atom in rdmol.GetAtoms()
                if is_port(atom)
    ]

neighbor_ports     = atomwise._neighbor_factory_by_condition(    condition=is_port)
has_neighbor_ports = atomwise._has_neighbor_factory_by_condition(condition=is_port)

# PORT MODIFICATION
@optional_in_place
def hydrogenate_rdmol_ports(rdmol : RDMol) -> None:
    '''Replace all port atoms with hydrogens'''
    for port_id in get_port_ids(rdmol):
        rdmol.GetAtomWithIdx(port_id).SetAtomicNum(1)