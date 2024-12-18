'''Unit tests for DISCERNMENT-related functionality'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import pytest
from polymerist.genutils.importutils.pkginspect import get_file_path_within_package
from polymerist.tests import data as testdata

import json
from pathlib import Path

from polymerist.genutils.sequences.discernment.inventory import SymbolInventory
from polymerist.genutils.sequences.discernment.strategies import DISCERNMENTStrategy


# DEFINE/LOAD HARD-CODED INPUTS AND EXPECTED OUTPUTS TO A PARTICULAR DISCERNMENT PROBLEM
@pytest.fixture(scope='module')
def word() -> str:
    return 'accg'

@pytest.fixture(scope='module')
def choice_bins() -> str:
    return ('bbc','aced','bad','daea','fccce','g','abcd','fggegc')

@pytest.fixture(scope='module')
def symbol_inventory(choice_bins) -> SymbolInventory:
    return SymbolInventory.from_bins(choice_bins)

@pytest.fixture(scope='module')
def solution_path() -> Path:
    return get_file_path_within_package('correct_discernment_solution.json', testdata)

@pytest.fixture(scope='module') 
def correct_solution(solution_path) -> set[tuple[int, ...]]:
    with solution_path.open('r') as file:
        solution = set(
            tuple(indices)
                for indices in json.load(file)
        )
    return solution


@pytest.mark.parametrize("ignore_multiplicities,unique_bins,DSClass", [(False, False, DSClass) for DSClass in DISCERNMENTStrategy.__subclasses__()]) # TODO: produce solutions with unique bins AND ignored multiplicities to fully test
class TestDISCERNMENTStrategies:
    all_results : dict[str, set[int]] = {} # cache solutions to avoid tedoius recalculations
    def test_preserves_symbol_inventory(
        self,
        word : str,
        symbol_inventory : SymbolInventory,
        correct_solution : set[tuple[int, ...]],
        ignore_multiplicities : bool,
        unique_bins : bool,
        DSClass : type[DISCERNMENTStrategy],
    ) -> None:
        '''Check to ensure that all implementations of DISCERNMENT solution strageties yield the same outputs and don't modify a provided symbol inventory
        Raises failure-specific Exception if inconsistency is detected, terminates silently (no Exception, returns None) otherwise'''
        mod_sym_inv = symbol_inventory.deepcopy() # create a copy of the symbol inventory to ensure any errant modifications do not affect other tests

        method_name = DSClass.__name__
        ds_strat = DSClass()
        proposed_solution = set(
            idxs
                for idxs in ds_strat.enumerate_choice_labels(
                    word,
                    mod_sym_inv,
                    ignore_multiplicities=ignore_multiplicities,
                    unique_bins=unique_bins
                )
        )
        self.all_results[method_name] = proposed_solution # cache for comparison in later tests
        assert (mod_sym_inv == symbol_inventory), f'Algorithm {method_name} does not produce to correct enumeration of bin labels'

    def test_solution_is_correct(
        self,
        word : str,
        symbol_inventory : SymbolInventory,
        correct_solution : set[tuple[int, ...]],
        ignore_multiplicities : bool,
        unique_bins : bool,
        DSClass : type[DISCERNMENTStrategy],
    ) -> None:
        method_name = DSClass.__name__
        assert (self.all_results[method_name] == correct_solution), f'Algorithm {method_name} does not return symbol inventory to original state after completion'

    def test_solution_strategies_are_consistent(
        self,
        word : str,
        symbol_inventory : SymbolInventory,
        correct_solution : set[tuple[int, ...]],
        ignore_multiplicities : bool,
        unique_bins : bool,
        DSClass : type[DISCERNMENTStrategy],
    ) -> None:
        method_name = DSClass.__name__
        proposed_solution = self.all_results[method_name]

        for other_method_name, other_solution in self.all_results.items(): # check against all other methods PRIOR TO INSERTION (minimal number of checks guaranteed to verify all pairwise checks)
            # check both symmetric differences to make sure no solution sequences are unique to either method 
            assert (proposed_solution - other_solution == set()) and (other_solution - proposed_solution == set()), f'Algorithms {method_name} and {other_method_name} produce inconsistent solutions'
