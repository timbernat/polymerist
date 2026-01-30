'''For generating linear polymer structure from monomer, sequence, and chain length information'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import logging
LOGGER = logging.getLogger(__name__)

import warnings
with warnings.catch_warnings(record=True): # suppress numerous and irritating mbuild deprecation warnings
    warnings.filterwarnings('ignore',  category=DeprecationWarning)
    from mbuild.lib.recipes.polymer import Polymer as MBPolymer

from typing import Optional

from .mbconvert import mbmol_from_mono_rdmol
from .sequencing import LinearCopolymerSequencer
from ..exceptions import MorphologyError
from ..monomers.repr import MonomerGroup
from ..estimation import estimate_n_atoms_linear
from ...genutils.textual.substrings import unique_string


def build_linear_polymer(
    monomers : MonomerGroup,
    n_monomers : int,
    sequence : str='A',
    minimize_sequence : bool=True,
    allow_partial_sequences : bool=False,
    add_Hs : bool=False,
    energy_minimize : bool=False,
    sequence_map : Optional[dict[str, str]]=None,
) -> MBPolymer:
    '''
    Builds a linear polymer structure from a specified pool of monomers, sequence, target chain length, and other parameters
    
    Parameters
    ----------
    monomers : MonomerGroup
        A group of fragments containing AT LEAST the distinct repeat units which occur in the target polymer
        
        IMPORTANT: if the "term_orient" field of the MonomerGroup is not set (with "head" and "tail" monomer designations),
        the first two terminal (1-valent) monomers in the group will be auto-assigned and taken as the head and tail, respectively,
        or, if there is only one terminal monomer present, it will be used as both the head and tail.
    n_monomers : int
        The number of monomer units in the target polymer chain
        
        This INCLUDES the terminal monomers in the count;
        E.g. n_monomers=10 with a head and tail group specified will induce 8 middle monomers
    sequence : str, default='A'
        A string of characters representing the sequence of MIDDLE monomers as they should occur within the polymer chain
        If the sequence is shorter than n_monomers - # end groups, the sequence will be repeated until the target chain length is reached.

        Each unique character in the string will be associated with a unique repeat unit
        in the provided MonomerGroup, as specified by `sequence_map` (see below)
        
        IMPORTANT: terminal monomers are not given by this sequence string, but are specified by either
        the "term_orient" field of the MonomerGroup or the auto-determined end groups if that is unset
    minimize_sequence : bool, default=True
        Whether to attempt to reduce the sequence provided into a minimal, repeating subsequence ("kernel")
        E.g. "ABABAB" will be reduced to 3*"AB" if this is set to True
        
        N.B.: this has NOTHING TO DO WITH energy minimization; that is controlled by the energy_minimize flag
    allow_partial_sequences : bool, default=False
        Whether to allow fractional repetitions of the sequence kernel to fill the target chain length
        
        For example, given a monomer group with head/tail specified and parameters
        n_monomers=10 and sequence="BAC" (inducing 10 - 2 = 8 middle monomers):
        * allow_partial_sequences=True will repeat the sequence 2 + 2/3 times,
          yielding the equivalent middle monomer sequence "BACBACBA", while
        * allow_partial_sequences=False would raise PartialBlockSequence Exception,
          since the sequence "BAC" cannot be repeated to fill 8 middle monomers exactly.
    add_Hs : bool, default=False
        Whether to instruct the mbuild Polymer recipe to cap uncapped terminal groups with hydrogens,
        in cases where the user has failed to provide ANY terminal monomers in the MonomerGroup
    energy_minimize : bool, default=False
        Whether to perform a brief UFF energy minimization after build to relax the resulting polymer structure
        Tends to give less-unphysical conformers for larger polymers but is significantly slower, especially for longer chains
    sequence_map : dict[str, str], default None
        Mapping from individual symbols (characters) in the sequence to the
        names of repeat units in monomers to associate with those symbols
        
        If no explicit mapping is provided, each unique symbol will be mapped to the
        MIDDLE repeat unit at that symbol's position when lexicographically sorted
        E.g. "BACA" will take the second, first, third, and first monomers defined
        in the MonomerGroup; so will "baca", "2131", "caea", and "tipi"
        
        For a more detailed example, suppose the contents of monomers.monomers looked like:
        >>> {
        >>>     'EG_term' : [...]
        >>>     'EG' : [...],
        >>>     'GA' : [...],
        >>>     'LA' : [...],
        >>>     ...
        >>> }
        
        Given sequence='bac' with sequence_map unset, a map of 
        >>> sequence_map = { # this is the default!!
        >>>     'a' : 'EG',
        >>>     'b' : 'GA',
        >>>     'c' : 'LA',
        >>> }
        would be generated by default (skipping over the terminal repeat unit 'EG_term') and
        the resulting polymer would be build as [head group]-[GA]-[EG]-[LA]-[GA]-[EG]-...
        
        On the other hand, supplying a map as:
        >>> sequence_map = {
        >>>     'a' : 'GA',
        >>>     'b' : 'LA',
        >>>     'c' : 'EG',
        >>> }
        would instead assemble the polymer as [head group]-[LA]-[GA]-[EG]-[LA]-[GA]-...
        
    Returns
    -------
    chain : mbuild.Compound
        An mbuild Polymer (Compound) object representing the assembled linear polymer chain
    '''
    # 1) DETERMINE NUMBER OF SEQUENCE REPEATS NEEDED TO MEET TARGET NUMBER OF MONOMER UNITS (IF POSSIBLE) - DEV: consider making a separate function
    end_groups = monomers.linear_end_groups() # cache end groups so they dont need to be recalculated when registering end groups
    end_group_names = [resname for (resname, _) in end_groups.values()]
    
    sequencer = LinearCopolymerSequencer(
        sequence_kernel=sequence,
        n_repeat_units=n_monomers,
        n_repeat_units_terminal=len(end_groups)
    )
    if minimize_sequence:
        sequencer.reduce() # identify minimal subsequences
    
    sequence_compliant, n_seq_repeats = sequencer.procrustean_alignment(allow_partial_sequences=allow_partial_sequences)
    LOGGER.info(
        f'Target chain length achievable with {sequencer.describe_tally()}, ' \
        f'namely with the sequence {sequencer.describe_order(end_group_names=end_group_names)}'
    )
    sequence_unique = unique_string(sequence_compliant, preserve_order=True) # only register a new monomer for each appearance of a new, unique symbol in the sequence
    
    # 2) REGISTER MONOMERS TO BE USED FOR CHAIN ASSEMBLY
    polymer = MBPolymer() 
    monomers_selected = MonomerGroup() # need only the subset of repeat units used in assembly to accurately assess linearity and estimate chain size (in current implementations)
    
    ## 2A) ADD MIDDLE MONOMERS TO CHAIN
    if sequence_map is None:
        sequence_map = {
            symbol : resname
                for symbol, (resname, middle_monomer) in zip(
                    sorted(sequence_unique),
                    monomers.iter_rdmols(term_only=False), # iterate over non-terminal repeat units only
                )
        }

    middle_monomer_table = monomers.rdmols(term_only=False) # cache this locally - DEV: ugly, but don't want to shift MonomerGroup API around too much
    if (num_symbols_unique := len(sequence_unique)) > (num_middle_monomers := len(middle_monomer_table)):
        raise ValueError(f'Too few unique repeat units ({num_middle_monomers}) to ascribe to each symbols of a {num_symbols_unique}-symbol sequence')
    
    for symbol in sorted(sequence_unique): # avoiding sequence_map.items() - won't count on sequence map being sorted (e.g. if a user passes one in)
        resname = sequence_map[symbol]     # another reason not to use sequence_map.items() here is we want a big, fat KeyError for undefined symbols in the sequence
        middle_monomer_choices = middle_monomer_table[resname]
        num_middle_monomer_choices : int = len(middle_monomer_choices)
        
        if num_middle_monomer_choices == 0:
            raise IndexError(f'No monomer templates for "{resname}" defined in MonomerGroup')
        elif num_middle_monomer_choices > 1:
            raise IndexError(f'Ambiguous choice for template "{resname}" ({num_middle_monomer_choices} templates defined for that label)')
        else:
            middle_monomer = middle_monomer_table[resname][0]
        
        LOGGER.info(f'Registering middle monomer {resname} (block identifier "{symbol}")')
        mb_monomer, linker_ids = mbmol_from_mono_rdmol(middle_monomer, resname=resname)
        polymer.add_monomer(compound=mb_monomer, indices=linker_ids)
        monomers_selected.monomers[resname] = monomers.monomers[resname]

    ## 2B) ADD TERMINAL MONOMERS TO CHAIN
    for head_or_tail, (resname, term_monomer) in end_groups.items():
        LOGGER.info(f'Registering terminal monomer {resname} (orientation "{head_or_tail}")')
        mb_monomer, linker_ids = mbmol_from_mono_rdmol(term_monomer, resname=resname)
        polymer.add_end_groups(compound=mb_monomer, index=linker_ids.pop(), label=head_or_tail, duplicate=False) # use single linker ID and provided head-tail orientation
        monomers_selected.monomers[resname] = monomers.monomers[resname]

    # 3) ASSEMBLE AND RETURN CHAIN
    if not monomers_selected.is_linear: # verify the selected monomers actually define a linear polymer
        raise MorphologyError('Linear polymer building does not support non-linear monomer input')
    
    n_atoms_est = estimate_n_atoms_linear(monomers_selected, n_monomers) # TODO: create new MonomerGroup with ONLY the registered monomers to guarantee accuracy
    LOGGER.info(f'Assembling linear {n_monomers}-mer chain (estimated {n_atoms_est} atoms)')
    
    polymer.build(n_seq_repeats, sequence=sequence_compliant, add_hydrogens=add_Hs) # "-2" is to account for term groups (in mbuild, "n" is the number of times to replicate just the middle monomers)
    LOGGER.info(f'Successfully assembled linear {n_monomers}-mer chain (exactly {polymer.n_particles} atoms)')
    
    # 4) OPTIONALLY, PERFORM FINAL UFF ENERGY MINIMIZATION
    if energy_minimize:
        LOGGER.info('Energy-minimizing chain to find more stable conformer')
        polymer.energy_minimize()
        LOGGER.info('Energy minimization completed')

    return polymer