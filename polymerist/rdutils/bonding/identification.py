'''Tools for determining how many and which bondable ports are in an RDKit Mol'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Generator, Iterable, Optional
from itertools import (
    combinations,
    product as cartesian_product,
)
from rdkit.Chem.rdchem import Mol

from .portlib import Port, get_ports, MolPortError
from ...genutils.iteration import iter_len


def get_total_port_degree_on_atom(rdmol : Mol, atom_id : int, target_flavor : Optional[int]=None):
    '''Get the maximum combined bond degree for all ports of a particular flavor on an atom (or of all flavors, if None is provided)'''
    return sum(
        int(port.bond.GetBondTypeAsDouble()) # sum over all matching-flavor ports on an atoms, weighting by bond order
            for port in get_ports(rdmol, target_atom_id=atom_id, target_flavor=target_flavor)
    )

def max_bondable_order_between_atoms(rdmol : Mol, atom_id_1 : int, atom_id_2 : int, target_flavor : int) -> int: # TODO : update or possible deprecate this
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

# def _get_bondable_ports_internal(rdmol : Mol) -> Generator[tuple[Port, Port], None, None]:
#     '''Generate all possible pairs of Ports within a single molecule which could be bonded together'''
#     return get_bondable_ports(rdmol, rdmol) # currently, this overcounts by exactly double, as a pair (p1, p2) and its copair (p2, p1) are considered distinct

def get_bondable_port_pairs(rdmol : Mol, atom_id_1 : Optional[int]=None, atom_id_2 : Optional[int]=None, flavor_pair : tuple[Optional[int], Optional[int]]=(None, None)) -> Generator[tuple[Port, Port], None, None]:
    '''Get every pair of ports within an RDKit Mol which could be bonded and match a specified flavor pair (could be partially or fully None for less specificity)
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

def get_first_bondable_port_pair(rdmol : Mol, atom_id_1 : Optional[int]=None, atom_id_2 : Optional[int]=None, flavor_pair : tuple[Optional[int], Optional[int]]=(None, None)) -> tuple[Port, Port]:
    '''Get the first pair of ports between atoms at given indices which are bondable, raising Exception if none exists'''
    try:
        return next(get_bondable_port_pairs(rdmol, atom_id_1=atom_id_1, atom_id_2=atom_id_2, flavor_pair=flavor_pair))
    except StopIteration:
        raise MolPortError(f'No bondable ports with flavors {flavor_pair} exist within Mol')

def get_num_bondable_port_pairs(rdmol : Mol, atom_id_1 : Optional[int]=None, atom_id_2 : Optional[int]=None, flavor_pair : tuple[Optional[int], Optional[int]]=(None, None)) -> int:
    '''Count how many bondable port pairs exist within a Mol'''
    return iter_len(get_bondable_port_pairs(rdmol, atom_id_1=atom_id_1, atom_id_2=atom_id_2, flavor_pair=flavor_pair))