'''Unit tests for DISCERNMENT-related functionality'''

import pytest
from polymerist.genutils.pkginspect import get_file_path_within_package
from polymerist.tests import data as testdata

import json
from pathlib import Path

from polymerist.genutils.sequences.discernment.inventory import SymbolInventory
from polymerist.genutils.sequences.discernment.strategies import DISCERNMENTStrategy


class DISCERNMENTInconsistencyError(Exception):
    '''Custom Exception for indicating inconsistencies between DISCERNMENT Strategy implementations'''
    pass


@pytest.fixture(scope='module')
def solution_path():
    return get_file_path_within_package('correct_discernment_solution.json', testdata)

@pytest.fixture(scope='module')
def correct_solution(solution_path):
    with solution_path.open('r') as file:
        solution = set(
            tuple(indices)
                for indices in json.load(file)
        )
    return solution


def test_discernment_algorithm_consistency(correct_solution : set[tuple[int, ...]], ignore_multiplicities : bool=False, unique_bins : bool=False) -> None:
    '''Check to ensure that all implementations of DISCERNMENT solution strageties yield the same outputs and don't modify a provided symbol inventory
    Raises failure-specific Exception if inconsistency is detected, terminates silently (no Exception, returns None) otherwise'''
    # hard-code inputs and expected solution
    WORD = 'accg'
    CHOICE_BINS = ('bbc','aced','bad','daea','fccce','g','abcd','fggegc')

    print(correct_solution)

    sym_inv = SymbolInventory.from_bins(CHOICE_BINS)
    orig_sym_inv = sym_inv.deepcopy()
    
    all_results : dict[str, set[int]] = {}
    for DSClass in DISCERNMENTStrategy.__subclasses__():
        method_name = DSClass.__name__
        ds_strat = DSClass()

        solution = set(
            idxs
                for idxs in ds_strat.enumerate_choice_labels(
                    WORD,
                    sym_inv,
                    ignore_multiplicities=ignore_multiplicities,
                    unique_bins=unique_bins
                )
        )
        if solution != correct_solution: # check that answer produces is accurate
            raise DISCERNMENTInconsistencyError(f'Algorithm {method_name} does not produce to correct enumeration of bin labels') 
        if sym_inv != orig_sym_inv: # check that the symbol inventory is unmodified
            raise DISCERNMENTInconsistencyError(f'Algorithm {method_name} does not return symbol inventory to original state after completion')
        for other_method_name, other_solution in all_results.items(): # check against all other methods PRIOR TO INSERTION (minimal number of checks guaranteed to verify all pairwise checks)
            # print(method_name, other_method_name)
            if (solution - other_solution != set()) or (other_solution - solution != set()): # check both symmetric differences to make sure none are 
                raise DISCERNMENTInconsistencyError(f'Algorithms {method_name} and {other_method_name} produce inconsistent solutions')

        # implicit else
        all_results[method_name] = solution