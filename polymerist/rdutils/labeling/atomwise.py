'''For obtaining info from and for labelling individual RDKit Atoms'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Callable, Generator
from rdkit.Chem.rdchem import Atom

# NEIGHBOR ATOM INFO
def get_num_bonds(atom : Atom) -> int:
    '''Returns number of explicit bonded connections an atom has (distinct from atom.GetExplicitValence due to bond orders)'''
    return len(atom.GetBonds())

def _get_neighbor_factory_by_condition(condition : Callable[[Atom], bool]) -> Callable[[Atom], Generator[Atom, None, None]]:
    '''Factory function for generating neighbor-search functions over Atoms by a boolean condition'''
    def neighbors_by_condition(atom : Atom) -> Generator[Atom, None, None]:
        '''Generate all neighboring atoms satisfying a condition'''
        for nb_atom in atom.GetNeighbors():
            if condition(nb_atom):
                yield nb_atom

    return neighbors_by_condition

def _has_neighbor_factory_by_condition(condition : Callable[[Atom], bool]) -> Callable[[Atom], bool]:
    '''Factory function for generating neighbor-search functions over Atoms by a boolean condition'''
    def has_neighbors_by_condition(atom : Atom) -> bool:
        '''Identify if any neighbors of an atom satisfy some condition'''
        return any(
            condition(nb_atom)
                for nb_atom in atom.GetNeighbors()
        )

    return has_neighbors_by_condition