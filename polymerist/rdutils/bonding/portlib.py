'''Specification for representing, defining, and characterizing selective intermolecular bond placeholders ("ports")'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import ClassVar, Generator, Optional
from dataclasses import dataclass, field

from rdkit import Chem
from rdkit.Chem.rdchem import Atom, Bond, Mol

from ...genutils.iteration import iter_len
from ...genutils.containers import UnorderedRegistry


# PORT DEFINITIONS
PORT_SCHEMA = '[#0D1]~[!#0]' # an atom of null type bonded (with exact one of any type of bond) to any non-null atom
PORT_QUERY  = Chem.MolFromSmarts(PORT_SCHEMA) # prototype molecule for a port

class MolPortError(Exception):
    '''Raised when port-related errors as encountered'''
    pass

def is_linker(rdatom : Atom) -> bool:
    '''Indicate whether an atom is a linker (intermonomer "*" type atom)'''
    return rdatom.GetAtomicNum() == 0

def get_num_linkers(rdmol : Mol) -> int:
    '''Count how many wild-type inter-molecule linker atoms are in a Mol'''
    return sum(
        is_linker(atom)
            for atom in rdmol.GetAtoms()
    )

@dataclass(frozen=True)
class Port:
    '''Class for encapsulating the components of a "port" bonding site (linker-bond-bridgehead)'''
    linker     : Atom
    bond       : Bond
    bridgehead : Atom

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

# PORT ENUMERATION
def get_port_ids(rdmol : Mol) -> Generator[tuple[int, int], None, None]:
    '''Get the linker and bridgehead indices of all ports found in an RDKit Mol'''
    for (linker_id, bh_id) in rdmol.GetSubstructMatches(PORT_QUERY):
        yield linker_id, bh_id # unpacked purely for self-documentation

def get_linker_ids(rdmol : Mol) -> Generator[int, None, None]:
    '''Get indices of all atoms which are ports'''
    for (linker_id, bh_id) in get_port_ids(rdmol):
        yield linker_id

def get_ports(rdmol : Mol, target_atom_id : Optional[int]=None, target_flavor : Optional[int]=None) -> Generator[Port, None, None]:
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

def get_num_ports(rdmol : Mol, target_atom_id : Optional[int]=None, target_flavor : Optional[int]=None) -> int: 
    '''Counts the number of port atoms present in a Mol'''
    return iter_len(get_ports(rdmol, target_atom_id=target_atom_id, target_flavor=target_flavor))

def get_single_port(rdmol : Mol) -> Port:
    '''Get the singular port of a Mol which contains only 1 port'''
    num_ports = get_num_ports(rdmol)
    if num_ports == 0:
        raise MolPortError("Can't return port for molecule which has no ports")
    
    if num_ports > 1:
        raise MolPortError("Molecule has multiple ports; choice of single port is ambiguous")
    
    return next(get_ports(rdmol)) # return if port count checks pass