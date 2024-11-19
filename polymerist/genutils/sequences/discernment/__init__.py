'''
Tools for solving the DISCERNMENT (Determination of Index Sequences from Complete Enumeration of Ransom Notes - Multiset Extension with Nonlexical Types) problem

DISCERNMENT problem definition:
    Given a "word" (a sequence of N symbols of type T), and a mapped sequence of "bins" (ordered collection of multisets of type T, each assigned a label of type L),
    enumerate all N-tuples of labels such that the symbols of the words could be drawn from the bins with those labels in that order
'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from .enumeration import DISCERNMENTSolver
from .strategies import (
    DISCERNMENTStrategyStack,
    DISCERNMENTStrategyCartesian,
    DISCERNMENTStrategyRecursive,
)