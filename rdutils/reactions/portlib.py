'''Utilities for finding, creating, and manipulating bonding sites (AKA "ports")'''

from typing import Generator
from rdkit import Chem

from ..rdtypes import RDAtom, RDMol
from ..labeling import atomwise

from ...genutils.decorators.functional import optional_in_place
from ...genutils.iteration import iter_len

# PORT INFO QUERY
PORT = Chem.AtomFromSmarts('[*]') # cast as lambda to generate new atom on each call (prevents cross-mutability issues)

def is_port(atom : RDAtom) -> bool:
    '''Shorthand for deciding if an atom is a bonding site ("port")''' 
    return atom.GetAtomicNum() == 0 # NOTE : can't check equality with PORT atom query (would match every atom)

def get_ports(rdmol : RDMol) -> Generator[RDAtom, None, None]:
    '''Generate all atoms which are ports'''
    for atom in rdmol.GetAtoms():
        if is_port(atom):  # NOTE : can't check equality with PORT atom query (would match every atom)
            yield atom

def get_port_ids(rdmol : RDMol) -> list[int]:
    '''Get indices of all atoms which are ports'''
    return [atom.GetIdx() for atom in get_ports(rdmol)]

def num_ports(rdmol : RDMol) -> int:
    '''Counts the number of port atoms present in a Mol'''
    return iter_len(get_ports(rdmol)) # NOTE : purposely avoid use of "len" to permit general iterators

neighbor_ports     = atomwise._neighbor_factory_by_condition(    condition=is_port)
has_neighbor_ports = atomwise._has_neighbor_factory_by_condition(condition=is_port)

# PORT MODIFICATION
@optional_in_place
def hydrogenate_rdmol_ports(rdmol : RDMol) -> None:
    '''Replace all port atoms with hydrogens'''
    for port_id in get_port_ids(rdmol):
        rdmol.GetAtomWithIdx(port_id).SetAtomicNum(1)