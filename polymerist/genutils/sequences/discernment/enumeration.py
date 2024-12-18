'''Front-facing solver facade for the DISCERNMENT Problem'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Generator, Iterable, Mapping, Sequence, Union

from .inventory import SymbolInventory, T, L
from .strategies import DISCERNMENTStrategy, DISCERNMENTStrategyStack


# TODO : add custom Exceptions to provide more detailed feedback?
class DISCERNMENTSolver():
    '''Encapsulation class for solving generalized ransom-note index enumeration problems for arbitrary words, bins, and solving algorithms'''
    def __init__(self, symbol_inventory : Union[SymbolInventory[T, L], Sequence[Iterable[T]], Mapping[L, Sequence[Iterable[T]]]], strategy : DISCERNMENTStrategy=DISCERNMENTStrategyStack()) -> None:
        if isinstance(symbol_inventory, SymbolInventory):
            self._symbol_inventory = symbol_inventory
        else:
            self._symbol_inventory = SymbolInventory.from_bins(symbol_inventory)
        self.strategy = strategy

    @property
    def symbol_inventory(self) -> SymbolInventory[T, L]:
        '''Cast symbol inventory as copy to avoid mutation during partial traversals (i.e. peek at first item)'''
        return self._symbol_inventory.deepcopy()

    def enumerate_choices(self, word : Sequence[T], ignore_multiplicities : bool=False, unique_bins : bool=False) -> Generator[L, None, None]:
        '''Enumerate all possible choices using specified solution strategy'''
        return self.strategy.enumerate_choice_labels(
            word=word,
            symbol_inventory=self.symbol_inventory,
            ignore_multiplicities=ignore_multiplicities,
            unique_bins=unique_bins,
        )

    def choice_solutions_exist(self, word : Sequence[T], ignore_multiplicities : bool=False, unique_bins : bool=False) -> bool:
        '''Precheck to see if a solution exists without attemting full enumeration'''
        if not self.symbol_inventory.contains_word(word, ignore_multiplicities=ignore_multiplicities):
            return False
        try:
            first_solution = next(self.enumerate_choices(word, ignore_multiplicities=ignore_multiplicities, unique_bins=unique_bins))
            return (first_solution != tuple()) # consider the null index sequence of the empty tuple to be a failure
        except StopIteration:
            return False
    choices_exist = choice_exists = solutions_exist = solution_exists = choice_solutions_exist # aliases for convenience
        
