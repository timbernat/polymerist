'''For obtaining info from and for labelling individual RDKit Atoms'''

from typing import Callable, Generator
from ..rdtypes import RDAtom


# NEIGHBOR ATOM INFO
def get_num_bonds(atom : RDAtom) -> int:
    '''Returns number of explicit bonded connections an atom has (distinct from atom.GetExplicitValence due to bond orders)'''
    return len(atom.GetBonds())

def _get_neighbor_factory_by_condition(condition : Callable[[RDAtom], bool]) -> Callable[[RDAtom], Generator[RDAtom, None, None]]:
    '''Factory function for generating neighbor-search functions over RDAtoms by a boolean condition'''
    def neighbors_by_condition(atom : RDAtom) -> Generator[RDAtom, None, None]:
        '''Generate all neighboring atoms satisfying a condition'''
        for nb_atom in atom.GetNeighbors():
            if condition(nb_atom):
                yield nb_atom

    return neighbors_by_condition

def _has_neighbor_factory_by_condition(condition : Callable[[RDAtom], bool]) -> Callable[[RDAtom], bool]:
    '''Factory function for generating neighbor-search functions over RDAtoms by a boolean condition'''
    def has_neighbors_by_condition(atom : RDAtom) -> bool:
        '''Identify if any neighbors of an atom satisfy some condition'''
        return any(
            condition(nb_atom)
                for nb_atom in atom.GetNeighbors()
        )

    return has_neighbors_by_condition