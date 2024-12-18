'''For calculating the edit distance between sequences and inspecting the edits needed to go between them'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Any, Generator, Sequence, Type, TypeVar, TypeAlias
from dataclasses import dataclass, field, replace
T = TypeVar('T')

from enum import Enum
import numpy as np

from ...typetools.numpytypes import N, M, Shape
from ...bits import int_to_bits


# REPRESENTATION CLASSES
class EditOperation(Enum): # TODO: reimplement as bitwise flags
    '''For annotating distinct kinds of sequence edits and their associated index offsets'''
    NULL         = 0 # NOTE: order of fields here intentional and CANNOT be modified!
    INSERTION    = 1 # ... this is because the bits of each index are the row and column
    DELETION     = 2 # ... offsets corresponding to that edit operation in a Wagner-Fischer matrix
    SUBSTITUTION = 3

    @property
    def bits(self) -> tuple[int, int]:
        '''Convert the integer value of the Enum field into its binary bits'''
        return tuple(int_to_bits(self.value, num_bits=2, as_list=True))
    offsets = bits # alias to make dependent code more readable

@dataclass
class EditInfo:
    '''for bundling together information about a sequence edit step'''
    edit_op  : EditOperation
    indices  : tuple[int, int]
    distance : int


# WAGNER-FISCHER MATRIX OPERATION
def compute_wf_matrix(seq1 : Sequence[T], seq2 : Sequence[T], int_type : Type=int) -> np.ndarray[Shape[N, M], int]:
    '''Compute (N+1)x(M+1) matrix of Levenshtein distances between all partial prefices of a pair of sequences
    where N and M are the lengths of the first and second sequence, respectively. Implements the Wagner-Fischer algorithm'''
    n, m = len(seq1), len(seq2)
    n_aug, m_aug = n + 1, m + 1

    # initialize matrix with all zeros apart from first row and column,
    wf_matrix = np.zeros((n_aug, m_aug), dtype=int_type)
    wf_matrix[:, 0] = np.arange(n_aug, dtype=int_type) # index along first column is same as number of deletions to get to empty sequence - NOTE: element [0, 0] overlaps here
    wf_matrix[0, :] = np.arange(m_aug, dtype=int_type) # index along first row is same as number of insertions to get from empty sequence - NOTE: element [0, 0] overlaps here

    # populate matrix by iterating on distances between sub-sequence problems
    for i, elem1 in enumerate(seq1, start=1):
        for j, elem2 in enumerate(seq2, start=1):
            wf_matrix[i, j] = 1 + np.min([
                wf_matrix[i-1, j],  # deletion at end of second sequence (after augmentation to end of first sequence)
                wf_matrix[i, j-1],  # insertion at end of second sequence
                wf_matrix[i-1, j-1] - int(elem1 == elem2) # substition of last elements between sequences (if elements are equal, then substitution costs nothing)
            ])
     # TODO: implement support for transposition weighting (i.e. Damerau-Levenshtein distance)

    return wf_matrix

def traverse_wf_matrix(wf_matrix : np.ndarray[Shape[N, M], int], begin_idxs : tuple[int, int]=(0, 0), end_idxs : tuple[int, int]=(-1, -1)) -> Generator[list[EditInfo], None, None]:
    '''Takes a Wagner-Fischer Levenshtein distance matrix and returns the indices of the minimal path through the matrix
    from the origin (i.e. empty sequences) to the '''
    assert(wf_matrix.ndim == 2)
    n_aug, m_aug = wf_matrix.shape

    if end_idxs == begin_idxs:
        yield [ # need to wrap in list for consistent typing in recursive calls
            EditInfo( # base case
                edit_op=EditOperation.NULL,
                indices=begin_idxs,
                distance=0
            )
        ]
    else:
        i, j = end_idxs
        i %= n_aug # ensure values are positive
        j %= m_aug # ensure values are positive
        curr_dist = wf_matrix[i, j]

        prev_edits = []
        for edit_op in EditOperation:
            if edit_op == EditOperation.NULL:
                continue # need to skip over to avoid RecursionError - TOSELF: would be nice to reimplement EditOperations as bitwise flags to streamline this

            di, dj = edit_op.offsets # unpack index offsets
            i_prev, j_prev = i - di, j - dj
            prev_edits.append(
                EditInfo(
                    edit_op=edit_op,
                    indices=(i_prev, j_prev),
                    distance=wf_matrix[i_prev, j_prev]
                )
            )

        min_prev_dist = min(ei.distance for ei in prev_edits)
        for edit_info in prev_edits:
            ret_edit_info = replace(edit_info, indices=(i, j)) # create a copy of the current edit info which reports the current 
            if (edit_info.edit_op == EditOperation.SUBSTITUTION) and (curr_dist == min_prev_dist):
                ret_edit_info.edit_op = EditOperation.NULL
                
            if edit_info.distance == min_prev_dist:
                for edit_steps in traverse_wf_matrix(wf_matrix, begin_idxs=begin_idxs, end_idxs=edit_info.indices): # recursive tail call through all possible predecessors
                    yield edit_steps + [ret_edit_info]

def describe_edits(seq1 : Sequence[T], seq2 : Sequence[T], int_type : Type=int, indicator : str=' -> ', delimiter : str='\n') -> Generator[str, None, None]:
    '''Describes step-by-step the insertion, deletion, or substitution operations needed to transform one sequence into another'''
    seqs = (seq1, seq2) # need to bundle for zipping later
    wf_matrix = compute_wf_matrix(seq1, seq2, int_type=int_type)

    for edits in traverse_wf_matrix(wf_matrix):
        edit_descs : list[str] = []
        for edit_info in edits:
            edit_op = edit_info.edit_op
            if edit_op == EditOperation.NULL:
                continue

            elem_edit_str = indicator.join( # TODO : add printout for how the word being modified looks after each change (notnjust which symbols are changing)
                str(seq[idx - 1] if to_show else None) # need to subtract 1 to get index of step previous to current edit
                    for idx, to_show, seq in zip(edit_info.indices, edit_op.bits, seqs)
            )

            edit_descs.append(f'{edit_info.edit_op.name}: {elem_edit_str}')
        yield delimiter.join(edit_descs)


# EDIT DISTANCES
def levenshtein_distance(seq1 : Sequence[T], seq2 : Sequence[T], int_type : Type=int) -> int:
    '''Compute the Levenshtein (edit) distance between a pair of sequences with elements of compatible type
    Denotes the minimal number of insertion, deletion, or substitution operations needed to transform either sequence into the other'''
    return compute_wf_matrix(seq1, seq2, int_type=int_type)[-1, -1]