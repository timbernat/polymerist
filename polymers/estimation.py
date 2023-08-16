'''For estimating properties of chains based on their constituent monomers and chain info'''

import numpy as np
from rdkit import Chem

from .exceptions import InsufficientChainLengthError
from ..genutils.iteration import iter_len
from ..monomers.repr import MonomerGroup
from ..rdutils.reactions.ports import num_ports as get_num_ports


def estimate_chain_len_linear(monomers : MonomerGroup, DOP : int) -> int:
    '''Given a set of monomers and the desired degree of polymerization, estimate the length of the resulting chain
    !NOTE! : As-implemented, only works for linear homopolymers and block copolymers with equal an distribution of monomers'''
    # TOSELF : omitted logging for now, as it gets repeated on EVERY cycle in when called estimate_DOP_lower
    num_mono = iter_len(monomers.iter_rdmols)

    mono_term    = np.zeros(num_mono, dtype=bool) # terminality of each monomer (i.e. whether or not it is a term group)
    mono_multip  = np.zeros(num_mono, dtype=int) # multiplicity of each polymer (i.e. how many times is occurs in a chain)
    mono_contrib = np.zeros(num_mono, dtype=int) # contribution of each monomer (i.e. how many atoms does it add to the chain)

    for i, (resname, monomer) in enumerate(monomers.iter_rdmols):
        num_atoms = monomer.GetNumAtoms()
        num_ports = get_num_ports(monomer)
        is_term = MonomerGroup.is_terminal(monomer)

        mono_term[i] = is_term
        mono_multip[i] = is_term # temporarily set middle monomer contribution to 0
        mono_contrib[i] = num_atoms - num_ports

    num_term = sum(mono_term)
    num_mid  = num_mono - num_term # assumed that all monomers are either terminal or not
    mono_multip[~mono_term] = (DOP - num_term) / num_mid # naive assumption that all middle monomers contribute rest of chain equally (for homopolymers, this is always true)

    N = mono_contrib @ mono_multip # compute dot product to yield final count
    
    return N

def estimate_DOP_lower(monomers : MonomerGroup, max_chain_len : int, min_DOP : int=3) -> int:
    '''Returns the largest DOP for a set of monomers which yields a chain no longer than the specified chain length'''
    base_chain_len = estimate_chain_len_linear(monomers, min_DOP)
    if base_chain_len > max_chain_len: # pre-check when optimization is impossible
        raise InsufficientChainLengthError(f'Even shortest possible chain (DOP={min_DOP}, N={base_chain_len}) is longer than the specified max length of {max_chain_len} atoms')

    DOP = min_DOP 
    while estimate_chain_len_linear(monomers, DOP + 1) < max_chain_len: # check if adding 1 more monomer keeps the length below the threshold
        DOP += 1

    return DOP

def estimate_DOP_upper(monomers : MonomerGroup, min_chain_len : int, min_DOP : int=3) -> int: # NOTE : as currently defined, this also subsumes the case when the estimate and calculated length are exactly equal
    '''Returns the smallest DOP for a set of monomers which yields a chain no shorter than the specified chain length'''
    return estimate_DOP_lower(monomers, min_chain_len, min_DOP=min_DOP) + 1 # by definition, this is just 1 monomer longer than the lower bound