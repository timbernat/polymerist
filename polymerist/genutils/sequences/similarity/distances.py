'''Implementations of calculation methods for sequence distance ("inverse similarity") metrics'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

# NOTE: much of this could be supplanted in the future by the well-implemented textdistance library (https://github.com/life4/textdistance)
from typing import Sequence, TypeVar
T = TypeVar('T')

from .edits import levenshtein_distance

# SEQUENCE METRICS
def hamming_distance(seq1 : Sequence[T], seq2 : Sequence[T]) -> int:
    '''Compute the Hamming distance between a pair of sequences with elements of compatible type (sequences must have the same length)
    Denotes the number of elements at the same positions in each sequence which are different'''
    if len(seq1) != len(seq2):
        raise ValueError('Cannot compute Hamming distance between sequences of different lengths')
    
    return sum(
        int(elem1 != elem2) # NOTE: type conversion not strictly necessary here, but done for self-documentation
            for elem1, elem2 in zip(seq1, seq2, strict=True)
    )

def jaccard_distance(seq1 : Sequence[T], seq2 : Sequence[T]) -> float:
    '''Compute the Jaccard distance between a pair of sequences with elements of compatible type
    Denotes the complement of the ratio of shared elements (intersection) to total elements (union)'''
    set1, set2 = set(seq1), set(seq2)
    size_intersect = len(set.intersection(set1, set2))
    size_union     = len(set.union(set1, set2))
    jaccard_coeff = size_intersect / size_union

    return 1 - jaccard_coeff
tanimoto_distance = jaccard_distance # TOSELF: debatable whether this alias is really accurate (literature suggests it may be context/field-dependent)

levenshtein_dist = edit_distance = edit_dist = levenshtein_distance # aliases for convenience