'''Tool for swapping bonds within and between RDKit Mols'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import logging
LOGGER = logging.getLogger(__name__)

from typing import Optional
from rdkit.Chem.rdchem import RWMol

from collections import Counter
from IPython.display import display # for Jupyter display support

from . import portlib
from .formation import increase_bond_order
from .dissolution import decrease_bond_order
from ..labeling import molwise

from ...genutils.sequences.seqops import int_complement
from ...genutils.decorators.functional import optional_in_place


def _is_valid_bond_derangement(bond_derangement : dict[int, tuple[int, int]]) -> bool:
    '''Determine whether an interatomic bond remapping describes a valid derangement'''
    # 1) check that each swap maps to a new element (i.e. no identity swaps)
    for begin_map_num, (curr_end_map_num, targ_end_map_num) in bond_derangement.items():
        if curr_end_map_num == targ_end_map_num:
            LOGGER.warn(f'Swap defined for initial index {begin_map_num} maps back to current partner ({curr_end_map_num} -> {targ_end_map_num})')
            return False
        
    # 2) check bijection (i.e. terminal atom remappings form a closed multiset)
    curr_end_counts, targ_end_counts = [Counter(i) for i in zip(*bond_derangement.values())] #  multisets are permissible for when multiple current/target bonds connect to the same atom 
    if curr_end_counts != targ_end_counts:
        LOGGER.warn('Bond derangement does not define a 1-1 correspondence between elements in the multiset')
        return False

    return True # only return if all above checks pass

@optional_in_place
def swap_bonds(rwmol : RWMol, bond_derangement : dict[int, tuple[int, int]], show_steps : bool=False) -> Optional[RWMol]:
    '''
    Takes a modifiable Mol and a bond derangement dict and performs the requested bond swaps
    Derangement dict should have th following form:
        keys   : int             = corresponds to the beginning atom of a bond
        values : tuple[int, int] = corresponds to the current end atom map number and target end atom map number (in that order) 
    Modifiable Mol can contain multiple disconnected molecular components
    ''' 
    # TODO : check for complete atom map num assignment
    if not _is_valid_bond_derangement(bond_derangement):
        raise ValueError('Invalid interatomic bond derangement provided')

    # determine non-interfering port flavors for new bonds (preserves parity between permutation sets)
    available_port_flavors = int_complement(molwise.get_isotopes(rwmol), bounded=False) # ensures newly-created temporary ports don't clash with any existing ones
    flavor_pair = (next(available_port_flavors), next(available_port_flavors)) # grab first two available flavors
    portlib.Port.bondable_flavors.insert(flavor_pair) # temporarily register pair as bondable

    # break current bonds
    for begin_map_num, (curr_end_map_num, _) in bond_derangement.items():
        decrease_bond_order(
            rwmol,
            *molwise.atom_ids_by_map_nums(rwmol, begin_map_num, curr_end_map_num),
            new_flavor_pair=flavor_pair,
            in_place=True # must be done in-place to allow optional_in_place decoration
        )

        if show_steps:
            print(f'{begin_map_num} --x-> {curr_end_map_num}')
            display(rwmol)

    # form new bonds - must be done AFTER breakage to ensure all necessary ports exist
    for begin_map_num, (_, targ_end_map_num) in bond_derangement.items():
        increase_bond_order(
            rwmol,
            *molwise.atom_ids_by_map_nums(rwmol, begin_map_num, targ_end_map_num), 
            flavor_pair=flavor_pair,
            in_place=True # must be done in-place to allow optional_in_place decoration
        )

        if show_steps:
            print(f'{begin_map_num} ----> {targ_end_map_num}')
            display(rwmol)

    # deregister bondable pair
    portlib.Port.bondable_flavors.pop(flavor_pair) 