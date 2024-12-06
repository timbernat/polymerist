'''For generating and manipulating sequences of symbols which correspond to monomer ordering in blocky and random copolymers'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import logging
LOGGER = logging.getLogger(__name__)

from ..exceptions import EndGroupDominatedChain, InsufficientChainLength, PartialBlockSequence
from ...genutils.textual.substrings import repeat_string_to_length


def procrustean_polymer_sequence_alignment(
        sequence : str,
        n_monomers_target : int,
        n_monomers_terminal : int,
        allow_partial_sequences : bool=False
    ) -> tuple[str, int]:
    '''
    For a given polymer block sequence "S", target linear chain length, and number of terminal monomers,
    Returns a sequence "P" and number of repeats "r" which, taken together, satisfy the following:
    - The number of monomers in r repeats of P plus the number of terminal monomers is precisely equal to the target number of monomers
    - The symbols in sequence P cycle through the symbols in S, in the order they appear in S
    - The number of times S is cycles through in P is always a rational multiple of the length of S
    If no satisfiable sequence-count pair can be found, raises an appropriate informative exception
    
    Named to reflect the fact that the original sequence S will be stretched or truncated to fit the given target sequence length
    
    Parameters
    ----------
    sequence : str
        A sequence indicating a periodic ordering of monomers in a linear polymer block (e.g. "A", "ABAC", etc)
        Each unique symbol in the sequence corresponds to a distinct monomer in the block
    n_monomers_target : int
        The desired number of monomers (including terminal monomers) in a polymer chain
    n_monomers_terminal : int
        The number of terminal monomers ("end groups") which are to be included in the chain
        in addition to the middle monomers described by "sequence"
    allow_partial_sequences : bool, default False
        Whether to allow fractional repeats of the original sequence in order to meet the target number of monomers
        
        For example, to construct a 12-mer chain with 2 end groups from the sequence "BACA", one would require 10 middle monomers
        which can only be achieved with 2.5 (10/4) sequence repeats, namely as "BACA|BACA|BA"; 

        This behavior may or may not be desired, depending on the use case, and can be controlled by this flag
    
    Returns
    -------
    sequence_procrustean : str
        A possibly modified version of the original polymer block sequence
    n_seq_repeats : int
        The number of times "sequence_procrustean" must be repeated to achieve the target sequence length
    
    Raises
    ------
    End GroupDominatedChain
        The number of terminal monomers exceed the number of total monomers
    PartialBlockSequence
        If a partial sequence repeat is required but disallowed (by setting allow_partial_sequences=False)
    InsufficientChainLength
        If the target number of monomers results in no middle monomers being included (i.e. neither full NOR partial sequence repeats)
    '''
    # Evaluate sizes of missing components from given values
    block_size = len(sequence)
    n_mono_middle = n_monomers_target - n_monomers_terminal # number of terminal monomers needed to reach target; in a linear chain, all monomers are either middle or terminal
    if n_mono_middle < 0:
        raise EndGroupDominatedChain(f'Registered number of terminal monomers exceeds requested chain length ({n_monomers_target}-mer chain can\'t possibly contain {n_monomers_terminal} terminal monomers)')
    
    n_seq_whole : int         # number of full sequence repeats to reach a number of monomers less than or equal to the target
    n_symbols_remaining : int # number of any remaining symbols in sequence (i.e. monomers) needed to close the gap to the target (allowed to be 0 if target is a multiple of the sequence length)
    n_seq_whole, n_symbols_remaining = divmod(n_mono_middle, block_size) 

    # Break down into cases by whether or not a whole number of sequence repeats is possible
    if n_symbols_remaining != 0: # a whole number of sequence repeats (including possibly 0) plus some fraction of a full block sequence
        if not allow_partial_sequences:
            raise PartialBlockSequence(
                f'Partial polymer block sequence required to meet target number of monomers ("{sequence[:n_symbols_remaining]}" prefix of sequence "{sequence}"). ' \
                'If this is acceptable, set "allow_partial_sequences=True" and try calling build routine again'
            )    
        sequence_procrustean = repeat_string_to_length(sequence, target_length=n_mono_middle, joiner='')
        n_seq_repeats = 1 # just repeat the entire mixed-fraction length sequence (no full sequence repeats to exploit)
    else: # for a purely-whole number of block sequence repeats
        if n_seq_whole < 1: # NOTE: if it were up to me, this would be < 0 to allow dimers, but mBuild has forced my hand
            raise InsufficientChainLength(
                f'{n_monomers_target}-monomer chain cannot accomodate both {n_monomers_terminal} end groups AND at least 1 middle monomer sequence'
            )
        sequence_procrustean = sequence # NOTE: rename here is for clarity, and for consistency with partial sequence case
        n_seq_repeats = n_seq_whole
        
    # Generate descriptive log message to summarize sequence modifications
    ## Determine info present for whole and partial sections
    desc_seq_counts_parts = []
    desc_seq_order_middle = []
    
    if n_seq_whole != 0: ## Whole sequence strings
        desc_seq_counts_parts.append(f'{n_seq_whole} whole {block_size}-sequence repeats')
        desc_seq_order_middle.append(f'{n_seq_whole}*[{sequence}]')
        
    if n_symbols_remaining != 0: ## Partial sequence strings
        desc_seq_counts_parts.append(f'a partial {n_symbols_remaining}/{block_size} sequence repeat')
        desc_seq_order_middle.append(f'[{sequence[:n_symbols_remaining]}]')
        
    ## Finalizing sequence counts descriptor parts
    tally_str = f'({n_seq_whole}*{block_size} + {n_symbols_remaining}) middle monomers + {n_monomers_terminal} terminal monomers = {n_monomers_target} total monomers)'
    if len(desc_seq_counts_parts) == 2:
        desc_seq_counts_parts.insert(1, ' and ') # include conjunction if a mixed (i.e. both whole and fractional) solution was found
    
    ## Finalizing sequence order descriptor parts
    desc_seq_order_parts = ['[END-GROUP]']*n_monomers_terminal # abut with correct amount of end group indicators
    desc_seq_order_parts[1:-1] = desc_seq_order_middle # insert middle sections for whole and partial sequences
    
    ## putting everything together
    LOGGER.info(f'Target chain length achievable with {"".join(desc_seq_counts_parts)};\n Namely, polymer will be sequenced as {" + ".join(desc_seq_order_parts)}, yielding {tally_str}')
        
    return sequence_procrustean, n_seq_repeats
