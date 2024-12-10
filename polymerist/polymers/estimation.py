'''For estimating properties of chains based on their constituent monomers and chain info'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import numpy as np

from .exceptions import InsufficientChainLength
from ..polymers.monomers.repr import MonomerGroup
from ..rdutils.bonding.portlib import get_num_ports


def estimate_n_atoms_linear(monomers : MonomerGroup, n_monomers : int) -> int:
    '''Given a set of monomers and the desired degree of polymerization, estimate the length of the resulting chain
    !NOTE! : As-implemented, only works for linear homopolymers and block copolymers with equal an distribution of monomers'''
    # TOSELF : omitted logging for now, as it gets repeated on EVERY cycle in when called estimate_n_monomers_supremum()
    num_mono = monomers.n_monomers
    mono_term    = np.zeros(num_mono, dtype=bool) # terminality of each monomer (i.e. whether or not it is a term group)
    mono_multip  = np.zeros(num_mono, dtype=int) # multiplicity of each polymer (i.e. how many times is occurs in a chain)
    mono_contrib = np.zeros(num_mono, dtype=int) # contribution of each monomer (i.e. how many atoms does it add to the chain)

    for i, (resname, monomer) in enumerate(monomers.iter_rdmols()):
        num_atoms = monomer.GetNumAtoms()
        num_ports = get_num_ports(monomer)
        is_term = MonomerGroup.is_terminal(monomer)

        mono_term[i] = is_term
        mono_multip[i] = is_term # temporarily set middle monomer contribution to 0
        mono_contrib[i] = num_atoms - num_ports

    num_term = sum(mono_term)
    num_mid  = num_mono - num_term # assumed that all monomers are either terminal or not
    mono_multip[~mono_term] = (n_monomers - num_term) / num_mid # naive assumption that all middle monomers contribute rest of chain equally (for homopolymers, this is always true)

    N = mono_contrib @ mono_multip # compute dot product to yield final count
    
    return N

def estimate_n_monomers_infimum(monomers : MonomerGroup, n_atoms_max : int, n_monomers_min : int=3) -> int:
    '''
    For a given collection of monomer fragments, returns the largest number of monomers which guarantees that
    a polymer chain made up of those monomers will have no more than the specified maximum number of atoms
    '''
    n_atoms_base = estimate_n_atoms_linear(monomers, n_monomers_min)
    if n_atoms_base > n_atoms_max: # pre-check when optimization is impossible
        raise InsufficientChainLength(f'Even shortest possible chain ({n_monomers_min} monomers, with {n_atoms_base} atoms) is longer than the specified max length of {n_atoms_max} atoms')

    n_monomers = n_monomers_min 
    while estimate_n_atoms_linear(monomers, n_monomers + 1) < n_atoms_max: # check if adding 1 more monomer keeps the length below the threshold
        n_monomers += 1

    return n_monomers

def estimate_n_monomers_supremum(monomers : MonomerGroup, n_atoms_min : int, n_monomers_min : int=3) -> int: # NOTE : as currently defined, this also subsumes the case when the estimate and calculated length are exactly equal
    '''
    For a given collection of monomer fragments, returns the smallest number of monomers which guarantees that
    a polymer chain made up of those monomers will have no fewer than the specified minimum number of atoms
    '''
    return estimate_n_monomers_infimum(monomers, n_atoms_min, n_monomers_min=n_monomers_min) + 1 # by definition, a ny more monomers than the infimum guarantees the chain will surpass a given number of atoms