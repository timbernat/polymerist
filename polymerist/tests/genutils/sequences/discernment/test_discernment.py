'''Unit tests for DISCERNMENT-related functionality'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import pytest

from typing import Type
from itertools import product as cartesian

from polymerist.genutils.importutils.pkginspect import get_file_path_within_package
from polymerist.tests import data as testdata

from polymerist.genutils.sequences.discernment.inventory import SymbolInventory
from polymerist.genutils.sequences.discernment.enumeration import DISCERNMENTSolver
from polymerist.genutils.sequences.discernment.strategies import (
    DISCERNMENTStrategy,
    DISCERNMENTStrategyCartesian,
    DISCERNMENTStrategyRecursive,
    DISCERNMENTStrategyStack,
)
from polymerist.genutils.sequences.discernment.examples import DISCERNMENTExample


# LOAD PREDEFINED INPUT-OUTPUT EXAMPLE PAIRS
## SHORT EXAMPLE - for lighter tests on testing aspects of solution besides correctness (i.e. side-effects)
EXAMPLE_SHORT = DISCERNMENTExample.from_file(
    get_file_path_within_package('DISCERNMENT/discernment_example_short.json', testdata)
)
params_short = [
    (strategy, *parameters)
        for strategy, parameters in cartesian(DISCERNMENTStrategy.__subclasses__(), EXAMPLE_SHORT.enumerate_test_inputs())
]

## LONG EXAMPLE - for more comprehensive tests on correctness
EXAMPLE_LONG  = DISCERNMENTExample.from_file(
    get_file_path_within_package('DISCERNMENT/discernment_example_long.json', testdata)
)
params_long = [
    (strategy, *parameters)
        for strategy, parameters in cartesian(DISCERNMENTStrategy.__subclasses__(), EXAMPLE_LONG.enumerate_test_inputs())
]


# TESTS PROPER
class TestDISCERNMENTStrategies:
    @pytest.mark.parametrize(
        'strategy_type,choice_bins,target_word,ignore_multiplicities,unique_bins,solution_expected',
        params_short, # no need to try lengthy examples here, just checking side-effects
    )
    def test_discernment_strategy_no_side_effects(
            self,
            strategy_type : Type[DISCERNMENTStrategy],
            choice_bins : tuple[list[str], ...],
            target_word : list[str],
            ignore_multiplicities : bool,
            unique_bins : bool,
            solution_expected : set[tuple[int, ...]]
        ) -> None:
        '''Test that the act of enumerating DISCERNMENT solutions leaves the symbol inventory undisturbed'''
        symbol_inventory = SymbolInventory.from_bins(choice_bins)
        symbol_inventory_ref = symbol_inventory.deepcopy() # create a copy of the symbol inventory to ensure any errant modifications do not affect other tests
        solver = DISCERNMENTSolver(symbol_inventory, strategy=strategy_type())
        
        for _ in solver.enumerate_choices(
                word=target_word,
                ignore_multiplicities=ignore_multiplicities,
                unique_bins=unique_bins,
            ):
                pass

        assert (symbol_inventory == symbol_inventory_ref), f'Algorithm {strategy_type.__name__} does not leave symbol inventory invariant'

    @pytest.mark.parametrize(
        'strategy_type,choice_bins,target_word,ignore_multiplicities,unique_bins,solution_expected',
        params_short + params_long, # try both long and short examples to maximize number of test cases
    )
    def test_discernment_strategy_soundness(
            self,
            strategy_type : Type[DISCERNMENTStrategy],
            choice_bins : tuple[list[str], ...],
            target_word : list[str],
            ignore_multiplicities : bool,
            unique_bins : bool,
            solution_expected : set[tuple[int, ...]]
        ) -> None:
        '''Test that DISCERNMENT solution strategy enumerates the correct label sequences for a given set of inputs'''
        symbol_inventory = SymbolInventory.from_bins(choice_bins)
        solver = DISCERNMENTSolver(symbol_inventory, strategy=strategy_type())
        
        solution_proposed = set(
            tuple(labels)
                for labels in solver.enumerate_choices(
                        word=target_word,
                        ignore_multiplicities=ignore_multiplicities,
                        unique_bins=unique_bins,
                    )
        )

        assert (solution_proposed == solution_expected), f'Algorithm {strategy_type.__name__} does not produce correct solutions to enumeration problem'
