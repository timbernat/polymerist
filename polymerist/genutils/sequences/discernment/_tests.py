'''Unit tests for DISCERNMENT-related functionality'''

# TODO: add logging

import json
from pathlib import Path

from .inventory import SymbolInventory
from .strategies import DISCERNMENTStrategy


class DISCERNMENTIncosistencyError(Exception):
    '''Custom Exception for indicating inconsistencies between DISCERNMENT Strategy implementations'''
    pass

def check_discernment_algorithm_consistency(ignore_multiplicities : bool=False, unique_bins : bool=False) -> None:
    '''Check to ensure that all implementations of DISCERNMENT solution strageties yield the same outputs and don't modify a provided symbol inventory
    Raises failure-specific Exception if inconsistency is detected, terminates silently (no Exception, returns None) otherwise'''
    # hard-code inputs and expected solution
    WORD = 'accg'
    CHOICE_BINS = ('bbc','aced','bad','daea','fccce','g','abcd','fggegc')

    LOCAL_DIR = Path(__file__).parent
    SOLUTION_PATH = LOCAL_DIR / 'correct_discernment_solution.json'

    with SOLUTION_PATH.open('r') as file:
        CORRECT_SOLUTION = set(
            tuple(indices)
                for indices in json.load(file)
        )

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
        if solution != CORRECT_SOLUTION: # check that answer produces is accurate
            raise DISCERNMENTIncosistencyError(f'Algorithm {method_name} does not produce to correct enumeration of bin labels') 
        if sym_inv != orig_sym_inv: # check that the symbol inventory is unmodified
            raise DISCERNMENTIncosistencyError(f'Algorithm {method_name} does not return symbol inventory to original state after completion')
        for other_method_name, other_solution in all_results.items(): # check against all other methods PRIOR TO INSERTION (minimal number of checks guaranteed to verify all pairwise checks)
            # print(method_name, other_method_name)
            if (solution - other_solution != set()) or (other_solution - solution != set()): # check both symmetric differences to make sure none are 
                raise DISCERNMENTIncosistencyError(f'Algorithms {method_name} and {other_method_name} produce inconsistent solutions')

        # implicit else
        all_results[method_name] = solution