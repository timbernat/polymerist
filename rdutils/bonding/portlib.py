'''Utilities for finding, creating, and manipulating bonding sites (AKA "ports")'''
from typing import ClassVar, Generator, Iterable, Optional

from itertools import combinations
from itertools import product as cartesian_product

from dataclasses import dataclass, field
from itertools import chain
from rdkit import Chem

from ..rdtypes import RDAtom, RDBond, RDMol
from ...genutils.iteration import iter_len
from ...genutils.containers import UnorderedRegistry
from ...genutils.decorators.functional import optional_in_place


# PORT DEFINITIONS
PORT_SCHEMA = '[#0D1]~[!#0]' # an atom of null type bonded (with exact one of any type of bond) to any non-null atom
PORT_QUERY  = Chem.MolFromSmarts(PORT_SCHEMA) # prototype molecule for a port

class MolPortError(Exception):
    '''Raised when port-related errors as encountered'''
    pass

def is_linker(rdatom : RDAtom) -> bool:
    '''Indicate whether an atom is a linker (intermonomer "*" type atom)'''
    return rdatom.GetAtomicNum() == 0

@dataclass(frozen=True)
class Port:
    '''Class for encapsulating the components of a "port" bonding site (linker-bond-bridgehead)'''
    linker     : RDAtom
    bond       : RDBond
    bridgehead : RDAtom

    bondable_flavors : ClassVar[UnorderedRegistry] = field(default=UnorderedRegistry((0, 0))) # by default, only two (0, 0)-flavor ("unlabelled") ports are bondable

    @property
    def flavor(self) -> int:
        '''Return the flavor of the port'''
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
            ((port_1.flavor, port_2.flavor) in Port.bondable_flavors)    # 1) port flavors are one of the registered pairs
            and port_1.bond.GetBondType() == port_2.bond.GetBondType()   # 2) the ports agree on the new bond order
            and port_1.bridgehead.GetIdx() != port_2.bridgehead.GetIdx() # 3) the two ports aren't connected to the same atom
            # and not port_1.bridgehead.Match(port_2.bridgehead) # TODO : workshop this comparison for query atoms TOSELF : this match is too general (matches by query, not identity)
        )
    
    def __hash__(self) -> int: # necessary to variation in seemingly identical RDKit atoms/bonds 
        '''Distills down unique info about a Port. Simplifies equality check and allows Ports to participate in sets as expected'''
        return hash((self.flavor, self.linker.GetIdx(), self.bridgehead.GetIdx(), self.bond.GetBondTypeAsDouble())) 
    
    def matches_flavor(self, target_flavor : Optional[int]) -> bool:
        '''Returns whether a port has a particular flavor'''
        if target_flavor is None:
            return True # always returning True if None is passed (simplifies non-specific defaults)
        
        return self.flavor == target_flavor


# PORT COUNTING AND INDEXING
# def get_num_ports(rdmol : RDMol) -> int: # NOTE : deprecated due to lack of selectivity for flavor
#     '''Counts the number of port atoms present in a Mol'''
#     return len(rdmol.GetSubstructMatches(PORT_QUERY))

def get_port_ids(rdmol : RDMol) -> Generator[tuple[int, int], None, None]:
    '''Get the linker and bridgehead indices of all ports found in an RDMol'''
    for (linker_id, bh_id) in rdmol.GetSubstructMatches(PORT_QUERY):
        yield linker_id, bh_id # unpacked purely for self-documentation

def get_linker_ids(rdmol : RDMol) -> Generator[int, None, None]:
    '''Get indices of all atoms which are ports'''
    for (linker_id, bh_id) in get_port_ids(rdmol):
        yield linker_id


# PORT ENUMERATION
def get_ports(rdmol : RDMol, target_atom_id : Optional[int]=None, target_flavor : Optional[int]=None) -> Generator[Port, None, None]:
    '''Find and generate all ports in a molecule. Can optioanlly narrow scope to port whose bridegehead is a particular atom and/or whose linker has a particular flavor'''
    target_atom = None if (target_atom_id is None) else rdmol.GetAtomWithIdx(target_atom_id)

    for (linker_id, bh_id) in get_port_ids(rdmol):
        port = Port(
            linker=rdmol.GetAtomWithIdx(linker_id),
            bond=rdmol.GetBondBetweenAtoms(linker_id, bh_id),
            bridgehead=rdmol.GetAtomWithIdx(bh_id)
        )

        if port.matches_flavor(target_flavor):
            # if (target_atom is None) or port.bridgehead.Match(target_atom):# will match if no atom is given OR if an atom is given and the bridgehead matches that Atom
            if (target_atom_id is None) or (port.bridgehead.GetIdx() == target_atom_id):# will match if no atom is given OR if an atom is given and the bridgehead matches that Atom
                yield port

def get_num_ports(rdmol : RDMol, target_atom_id : Optional[int]=None, target_flavor : Optional[int]=None) -> int: 
    '''Counts the number of port atoms present in a Mol'''
    return iter_len(get_ports(rdmol, target_atom_id=target_atom_id, target_flavor=target_flavor))

def get_single_port(rdmol : RDMol) -> Port:
    '''Get the singular port of a Mol which contains only 1 port'''
    num_ports = get_num_ports(rdmol)
    if num_ports == 0:
        raise MolPortError("Can't return port for molecule which has no ports")
    
    if num_ports > 1:
        raise MolPortError("Molecule has multiple ports; choice of single port is ambiguous")
    
    return next(get_ports(rdmol)) # return if port count checks pass

# PORT BOND-PAIR ENUMERATION
def get_total_port_degree_on_atom(rdmol : RDMol, atom_id : int, target_flavor : Optional[int]=None):
    '''Get the maximum combined bond degree for all ports of a particular flavor on an atom (or of all flavors, if None is provided)'''
    return sum(
        int(port.bond.GetBondTypeAsDouble()) # sum over all matching-flavor ports on an atoms, weighting by bond order
            for port in get_ports(rdmol, target_atom_id=atom_id, target_flavor=target_flavor)
    )

def max_bondable_order_between_atoms(rdmol : RDMol, atom_id_1 : int, atom_id_2 : int, target_flavor : int) -> int: # TODO : update or possible deprecate this
    '''Return the highest possible bond order which could be created between a pair of atoms''' # NOTE : can't just count number of bondable pairs, since these include all possible permutations
    if target_flavor is None:
        raise NotImplementedError # TODO : find better way to handle this (lack of selectivity with no flavor prevents accurate bondability estimation)
    
    return min( # limited by which atom has the fewest bonding sites
        get_total_port_degree_on_atom(rdmol, atom_id_1, target_flavor=target_flavor),
        get_total_port_degree_on_atom(rdmol, atom_id_2, target_flavor=target_flavor),
    )

def enumerate_bondable_port_pairs(ports_1 : Iterable[Port], ports_2 : Iterable[Port]) -> Generator[tuple[Port, Port], None, None]:
    '''Generate all possible pairs of Ports between two molecules which could be bonded together'''
    for port_1, port_2 in cartesian_product(ports_1, ports_2):
        if Port.are_bondable(port_1, port_2):
            yield (port_1, port_2)

def enumerate_bondable_port_pairs_internal(ports : Iterable[Port]) -> Generator[tuple[Port, Port], None, None]: # TODO : find way to merge this with enumerate_bondable_port_pairs
    '''Generate all possible pairs of Ports within a single molecule which could be bonded together'''
    for port_1, port_2 in combinations(ports, 2):
        if Port.are_bondable(port_1, port_2):
            yield (port_1, port_2)

# def _get_bondable_ports_internal(rdmol : RDMol) -> Generator[tuple[Port, Port], None, None]:
#     '''Generate all possible pairs of Ports within a single molecule which could be bonded together'''
#     return get_bondable_ports(rdmol, rdmol) # currently, this overcounts by exactly double, as a pair (p1, p2) and its copair (p2, p1) are considered distinct

def get_bondable_port_pairs(rdmol : RDMol, atom_id_1 : Optional[int]=None, atom_id_2 : Optional[int]=None, flavor_pair : tuple[Optional[int], Optional[int]]=(None, None)) -> Generator[tuple[Port, Port], None, None]:
    '''Get every pair of ports within an RDMol which could be bonded and match a specified flavor pair (could be partially or fully None for less specificity)
    Can optionally localize search to only have port bridegehead coinciding with atoms at particular indices
    Pairs are returned n descending order of linker atom index to simplify atom removal if modifying later'''
    # possible_port_pairs = enumerate_bondable_port_pairs_internal( # unpack values after iterating; TODO : make this internal to avoid double counting
    #     chain( # concatenate together two generators (one for each atom-flavor pair)
    #         *(get_ports(rdmol, target_atom_id=atom_id, target_flavor=flavor) # NOTE : need star unpacking here to pass ports onto thru chain (rather than just passing the two generators)
    #             for (atom_id, flavor) in zip((atom_id_1, atom_id_2), flavor_pair)
    #         )
    #     )
    # )
    possible_ports = set() # treated as set to account for uniqueness
    for atom_id in (atom_id_1, atom_id_2): # TODO : find less nested way to do this
        for flavor in flavor_pair:
            for port in get_ports(rdmol, target_atom_id=atom_id, target_flavor=flavor):
                possible_ports.add(port) 

    for (port_1, port_2) in enumerate_bondable_port_pairs_internal(possible_ports):
        if Port.are_bondable(port_1, port_2):
            if port_1.linker.GetIdx() >= port_2.linker.GetIdx(): # sort in-place to ensure highest-index linker is first (MATTER FOR ORDER OF ATOM REMOVAL)
                yield (port_1, port_2)
            else:
                yield (port_2, port_1)

def get_first_bondable_port_pair(rdmol : RDMol, atom_id_1 : Optional[int]=None, atom_id_2 : Optional[int]=None, flavor_pair : tuple[Optional[int], Optional[int]]=(None, None)) -> tuple[Port, Port]:
    '''Get the first pair of ports between atoms at given indices which are bondable, raising Exception if none exists'''
    try:
        return next(get_bondable_port_pairs(rdmol, atom_id_1=atom_id_1, atom_id_2=atom_id_2, flavor_pair=flavor_pair))
    except StopIteration:
        raise MolPortError(f'No bondable ports with flavors {flavor_pair} exist within Mol')

def get_num_bondable_port_pairs(rdmol : RDMol, atom_id_1 : Optional[int]=None, atom_id_2 : Optional[int]=None, flavor_pair : tuple[Optional[int], Optional[int]]=(None, None)) -> int:
    '''Count how many bondable port pairs exist within a Mol'''
    return iter_len(get_bondable_port_pairs(rdmol, atom_id_1=atom_id_1, atom_id_2=atom_id_2, flavor_pair=flavor_pair))
