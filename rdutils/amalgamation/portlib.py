'''Utilities for finding, creating, and manipulating bonding sites (AKA "ports")'''
from typing import Generator, Iterable

from itertools import combinations
from itertools import product as cartesian_product

from dataclasses import dataclass
from collections import defaultdict

from rdkit import Chem
from ..rdtypes import RDAtom, RDBond, RDMol
from ...genutils.decorators.functional import optional_in_place


# PORT DEFINITIONS
PORT_SCHEMA = '[#0D1]~[!#0]' # an atom of null type bonded (with exact one of any type of bond) to any non-null atom
PORT_QUERY  = Chem.MolFromSmarts(PORT_SCHEMA) # prototype molecule for a port

class MolPortError(Exception):
    '''Raised when port-related errors as encountered'''
    pass

@dataclass(frozen=True)
class Port:
    '''Class for encapsulating the components of a "port" bonding site (linker-bond-bridgehead)'''
    linker     : RDAtom
    bond       : RDBond
    bridgehead : RDAtom

    @property
    def desig(self) -> int:
        '''Return the designation of the port'''
        return self.linker.GetIsotope()
    
    def __eq__(self, other : 'Port') -> bool:
        if not isinstance(other, Port):
            raise TypeError(f'Cannot compare {self.__class__.__name__} to object of type {other.__class__.__name__}')
        
        return (
            self.linker.Match(other.linker) # must use explicit "Match" call to correctly compare QueryAtoms and QueryBonds (may be distinct but equivalent objects)
            and self.bond.Match(other.bond)
            and self.bridgehead.Match(other.bridgehead)
        )
    
    @staticmethod
    def are_bondable(port_1 : 'Port', port_2 : 'Port') -> bool:
        '''Determine if two port atoms can be combined into a new bond'''
        return (
            port_1.desig == port_2.desig                                 # 1) port are of the same designation (i.e. matching "flavor")
            and port_1.bond.GetBondType() == port_2.bond.GetBondType()   # 2) the ports agree on the new bond order
            and port_1.bridgehead.GetIdx() != port_2.bridgehead.GetIdx() # 3) the two ports aren't connected to the same atom
            # and not port_1.bridgehead.Match(port_2.bridgehead) # TODO : workshop this comparison for query atoms TOSELF : this match is too general (matches by query, not identity)
        )


# PORT COUNTING AND INDEXING
def get_num_ports(rdmol : RDMol) -> int:
    '''Counts the number of port atoms present in a Mol'''
    return len(rdmol.GetSubstructMatches(PORT_QUERY))

def get_port_ids(rdmol : RDMol) -> Generator[tuple[int, int], None, None]:
    '''Get the linker and bridgehead indices of all ports found in an RDMol'''
    for (linker_id, bh_id) in rdmol.GetSubstructMatches(PORT_QUERY):
        yield linker_id, bh_id # unpacked purely for self-documentation

def get_linker_ids(rdmol : RDMol) -> Generator[int, None, None]:
    '''Get indices of all atoms which are ports'''
    for (linker_id, bh_id) in get_port_ids(rdmol):
        yield linker_id


# PORT ENUMERATION
def get_ports(rdmol : RDMol) -> Generator[Port, None, None]:
    '''Find and generate all ports in a molecule'''
    for (linker_id, bh_id) in get_port_ids(rdmol):
        yield Port(
            linker=rdmol.GetAtomWithIdx(linker_id),
            bond=rdmol.GetBondBetweenAtoms(linker_id, bh_id),
            bridgehead=rdmol.GetAtomWithIdx(bh_id)
        )

def get_single_port(rdmol : RDMol) -> Port:
    '''Get the singular port of a Mol which contains only 1 port'''
    num_ports = get_num_ports(rdmol)
    if num_ports == 0:
        raise MolPortError("Can't return port for molecule which has no ports")
    
    if num_ports > 1:
        raise MolPortError("Molecule has multiple ports; choice of single port is ambiguous")
    
    return next(get_ports(rdmol)) # return if port count checks pass

# def get_ports_by_desig(rdmol : RDMol) -> dict[int, list[Port]]:
#     '''Returns dict of all port info for a Mol, keyed by port designation'''
#     ports_dict = defaultdict(list)
#     for port in get_ports(rdmol):
#         ports_dict[port.desig].append(port)
    
#     return ports_dict

def get_ports_with_desig(rdmol : RDMol, target_desig : int=0) -> Generator[Port, None, None]:
    '''Generate all ports in a molecule with a particular port designation'''
    for port in get_ports(rdmol):
        if port.desig == target_desig:
            yield port

def get_ports_on_atom_at_idx(rdmol : RDMol, atom_id : int) -> Generator[Port, None, None]:
    '''Generate all Ports in an RDMol whose bridegehead atom is the atom at the specified atomic index'''
    targ_atom = rdmol.GetAtomWithIdx(atom_id)
    for port in get_ports(rdmol):
        if port.bridgehead.Match(targ_atom):
            yield port


# PORT BOND-PAIR ENUMERATION
def get_bondable_port_pairs(ports_1 : Iterable[Port], ports_2 : Iterable[Port]) -> Generator[tuple[Port, Port], None, None]:
    '''Generate all possible pairs of Ports between two molecules which could be bonded together'''
    for port_1, port_2 in cartesian_product(ports_1, ports_2):
        if Port.are_bondable(port_1, port_2):
            yield (port_1, port_2)

# def _get_bondable_ports_internal(rdmol : RDMol) -> Generator[tuple[Port, Port], None, None]:
#     '''Generate all possible pairs of Ports within a single molecule which could be bonded together'''
#     return get_bondable_ports(rdmol, rdmol) # currently, this overcounts by exactly double, as a pair (p1, p2) and its copair (p2, p1) are considered distinct

def get_bondable_port_pairs_internal(ports : Iterable[Port]) -> Generator[tuple[Port, Port], None, None]:
    '''Generate all possible pairs of Ports within a single molecule which could be bonded together'''
    for port_1, port_2 in combinations(ports, 2):
        if Port.are_bondable(port_1, port_2):
            yield (port_1, port_2)


@optional_in_place # temporarily placed here for backwards-compatibility reasons
def hydrogenate_rdmol_ports(rdmol : RDMol) -> None:
    '''Replace all port atoms with hydrogens'''
    for port_id in get_port_ids(rdmol):
        rdmol.GetAtomWithIdx(port_id).SetAtomicNum(1)
