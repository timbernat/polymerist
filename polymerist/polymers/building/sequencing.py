'''For generating and manipulating sequences of symbols which correspond to monomer ordering in blocky and random copolymers'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import logging
LOGGER = logging.getLogger(__name__)

from typing import Iterable, Optional
from dataclasses import dataclass, field, asdict

from ...genutils.textual.substrings import shortest_repeating_substring, repeat_string_to_length
from ...genutils.fileutils.jsonio.jsonify import make_jsonifiable
from ..exceptions import EndGroupDominatedChain, InsufficientChainLength, EmptyBlockSequence, PartialBlockSequence


@make_jsonifiable
@dataclass
class LinearCopolymerSequencer:
    '''
    For encapsulating information about the sequence of repeat units in a periodic, linear copolymer
    Also covers, as trivial special cases, homopolymers and alternating copolymers
    
    Parameters
    ----------
    sequence_kernel : str
        A sequence indicating a periodic ordering of monomers in a linear polymer block (e.g. "A", "ABAC", etc)
        Each unique symbol in the sequence corresponds to a distinct monomer in the block
    n_repeat_units : int
        The desired total number of monomers (including terminal monomers) in a polymer chain
    n_monomers_terminal : int
        The number of terminal monomers ("end groups") which are to be included in the chain
        in addition to the middle monomers described by "sequence"
        
    Raises
    ------
    EmpyBlockSequence
        The sequence provided is empty (can't be used to define nonzero-length chain)
    End GroupDominatedChain
        The number of terminal monomers exceed the number of total monomers
    '''
    sequence_kernel : str
    n_repeat_units : int
    n_repeat_units_terminal : int = 0
    
    # Attribute checks and modifications
    def __post_init__(self) -> None:
        if not self.sequence_kernel:
            raise EmptyBlockSequence('Must provide non-empty sequence kernel to yield a valid (co)polymer sequence')
    
        if self.n_repeat_units_middle < 0:
            raise EndGroupDominatedChain(
                f'Number of terminal monomers exceeds requested chain length; ({self.n_repeat_units}-mer ' \
                f'chain can\'t possibly contain {self.n_repeat_units_terminal} terminal monomers)'
            )
            
    def copy(self) -> 'LinearCopolymerSequencer':
        '''Returns another equivalent instance of the current sequence info more efficiently than a complete deepcopy'''
        return self.__class__(**asdict(self))
            
    def reduce(self) -> None:
        '''
        Determines if there is a shorter repeating subsequence making up the current sequence kernel
        If there is, adjusts the sequence kernel to that minimal sequence; does nothing otherwise
        
        Reduction is idempotent, and guarantees that the smallest possible kernel is used when sequencing
        '''
        minimal_subsequence = shortest_repeating_substring(self.sequence_kernel)
        kernel_period = self.block_size // len(minimal_subsequence) # account for any periodic shortening WITHIN the kernel
        
        if kernel_period == 1:
            LOGGER.info(f'Sequence kernel "{self.sequence_kernel}" is already fully reduced; no changes made')
            return
        else:
            LOGGER.info(
                f'Sequence kernel "{self.sequence_kernel}" can be further decomposed as {kernel_period}*"{minimal_subsequence}"; ' \
                f'Setting kernel to minimal subsequence "{minimal_subsequence}"'
            )
            self.sequence_kernel = minimal_subsequence
    
    def reduced(self) -> 'LinearCopolymerSequencer':
        '''Return a sequence-reduced version of the current sequence info'''
        clone = self.copy()
        clone.reduce()
        
        return clone
        
    # Properties derived from sequence kernel and target chain lengths
   
    @property
    def n_repeat_units_middle(self) -> int:
        '''Number of middle (i.e. non-terminal) repeat units'''
        return self.n_repeat_units - self.n_repeat_units_terminal

    # Whole sequence periods
    @property
    def block_size(self) -> int:
        '''Number of repeat units units in one whole iteration of the kernel block'''
        return len(self.sequence_kernel)
    period = block_size
    
    @property
    def n_full_periods(self) -> int:
        '''
        Largest number of complete repetitions of the sequence kernel which, when taken
        together, contain no more repeats units than the specified number of middle units
        '''
        return self.n_repeat_units_middle // self.block_size
    
    # Partial sequence residues
    @property
    def n_residual_repeat_units(self) -> int:
        '''
        Difference between number of middle repeat units and units which
        would occur in maximal full periods of the kernel
        
        By construction, is no greater than the block size and is
        identically zero exactly when a whole number of kernel repeats
        '''
        return self.n_repeat_units_middle % self.block_size
    n_residual_symbols = n_res = n_residual_repeat_units
    
    @property
    def has_residual(self) -> bool:
        '''
        Whether or not the target number of middle repeat units
        can be attained by a whole number of kernel repeats
        '''
        return bool(self.n_residual_repeat_units)
    
    @property
    def sequence_residual(self) -> str:
        '''Partial repeat of the kernel sequence needed to attain the speficied number of middle units'''
        return self.sequence_kernel[:self.n_residual_repeat_units]
    residual = sequence_residual
    
    ## PROCRUSTEAN sequence alignment
    def procrustean_alignment(self, allow_partial_sequences : bool=False) -> tuple[str, int]:
        '''
        PROCRUSTEAN: Periodic Recurrence Of Cyclic Repeat Unit Sequences, Truncated to an Exact and Arbitrary Number
        Stretches or truncates the sequence kernel to achieve a target sequence length, cycling through the kernel's period as many times as needed
        
        Algorithm produces a sequence string "P" and number of repeats "r" which, taken together, satisfy the following:
        - The number of units in r repeats of P plus the number of terminal monomers is precisely equal to the target number of monomers
        - The units in P cycle through the units in S, in the order they appear in S
        - The number of times S is cycled through in P is always a rational multiple of the length of S
        If no satisfiable sequence-count pair can be found, raises an appropriate informative exception
        '''
        if not self.has_residual: # the case where the target length happens to consist of a whole-number of repeats of the kernel
            if self.n_full_periods < 1: # NOTE: if it were up to me, this would be < 0 to allow dimers, but mBuild has forced my hand
                raise InsufficientChainLength(
                    f'{self.n_repeat_units}-monomer chain cannot accomodate both {self.n_repeat_units_terminal} end groups AND at least 1 middle monomer sequence'
                )
            sequence_procrustean = self.sequence_kernel
            n_seq_repeats = self.n_full_periods
        else:
            if not allow_partial_sequences:
                raise PartialBlockSequence(
                    f'Partial polymer block sequence required to meet target number of monomers ("{self.residual}" prefix of sequence "{self.sequence_kernel}");\n' \
                    'If this is acceptable, set "allow_partial_sequences=True" and try calling build routine again'
                )    
            sequence_procrustean = repeat_string_to_length(self.sequence_kernel, target_length=self.n_repeat_units_middle, joiner='')
            n_seq_repeats = 1 # just repeat the entire mixed-fraction length sequence (no full sequence repeats to exploit)
            
        return sequence_procrustean, n_seq_repeats
    
    def describe_order(self, end_group_names : Optional[Iterable[str]]=None, default_end_group_name : str='END-GROUP') -> str:
        '''Descriptive string presenting a condensed view of the order of repeat units in the final sequence'''
        # Assign names for end groups
        if end_group_names is None:
            end_group_names = [f'[{default_end_group_name}]']*self.n_repeat_units_terminal 
        else:
            end_group_names = [f'[{end_group_name}]' for end_group_name in end_group_names] # unpack into list and enforce correct number of names
        if (num_names_provided := len(end_group_names)) != self.n_repeat_units_terminal: # DEV: consider supporting filling in missing names with default in future
            raise IndexError(f'Defined sequence info with {self.n_repeat_units_terminal} end groups, but only provided names for {num_names_provided}')
        
        # Insert middle omnomer parts as necessary
        sequence_middle = []
        if self.n_full_periods != 0: ## Whole sequence strings
            sequence_middle.append(f'{self.n_full_periods}*[{self.sequence_kernel}]')
        if self.has_residual: ## Partial sequence strings
            sequence_middle.append(f'[{self.residual}]')
            
        # Abut with correct amount of end group indicators
        sequence_parts = end_group_names[:] 
        sequence_parts[1:-1] = sequence_middle
        
        return ' + '.join(sequence_parts)
    
    def describe_tally(self) -> str:
        '''Descriptive string indicating how all parts of the overall sequence contribute to the target number of repeat units'''
        desc_seq_counts_parts = []
        if self.n_full_periods != 0: ## Whole sequence strings
            desc_seq_counts_parts.append(f'{self.n_full_periods} whole {self.block_size}-sequence repeat(s)')
        if self.has_residual: ## Partial sequence strings
            desc_seq_counts_parts.append(f'a partial {self.n_residual_repeat_units}/{self.block_size} sequence repeat')
            
        return ' and '.join(desc_seq_counts_parts)
        