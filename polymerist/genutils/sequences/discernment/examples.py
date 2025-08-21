'''Encapsulation for example input-output pairs to a DISCERNMENT problem; intended to facilitate unit testing'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Any, ClassVar, Generator, Iterable, Sequence, TypeVar
T = TypeVar('T')
L = TypeVar('L')

from dataclasses import dataclass
from collections import defaultdict

import json
from pathlib import Path


@dataclass
class DISCERNMENTExample:
    '''For encapsulating pre-made DISCERNMENT example input-output pairs'''
    choice_bins : tuple[Sequence[T], ...]
    target_word : Sequence[T]
    solutions : dict[tuple[bool, bool], set[tuple[int, ...]]]
    
    INDENT : ClassVar[int] = 4

    def to_json(self) -> dict[str, Any]:
        '''Write contents to JSON-serializable dict'''
        solutions_nested = defaultdict(dict) 
        for (ignore_multiplicities, unique_bins), solution in self.solutions.items():
            # nest flag keys, since JSON is too weak to handle list keys 
            # convert to int since JSON bools are not read as Python bools
            solutions_nested[int(ignore_multiplicities)][int(unique_bins)] = [list(i) for i in solution]

        return {
            'choice_bins' : [list(choice_bin) for choice_bin in self.choice_bins],
            'target_word' : list(self.target_word),
            'solutions' : dict(solutions_nested),
        }

    def to_file(self, example_path : Path) -> None:
        '''Save this example to a JSON file'''
        if isinstance(example_path, str):
            example_path = Path(example_path)

        with example_path.open('w') as file:
            json.dump(self.to_json(), file, indent=self.INDENT)

    @classmethod
    def from_file(cls, example_path: Path) -> 'DISCERNMENTExample':
        '''Load a DISCERNMENT example from a JSON file'''
        if isinstance(example_path, str):
            example_path = Path(example_path)

        with example_path.open('r') as file:
            json_dict = json.load(file)

        return cls(
            choice_bins=tuple(json_dict['choice_bins']),
            target_word=json_dict['target_word'],
            solutions={
                (bool(int(ignore_multiplicities)), bool(int(unique_bins))) : set(tuple(i) for i in solution)
                    for ignore_multiplicities, subsolns in json_dict['solutions'].items()
                        for unique_bins, solution in subsolns.items()
            }
        )
        
    def enumerate_test_inputs(self) -> Generator[
            tuple[
                tuple[Sequence[T], ...],
                Sequence[T],
                bool,
                bool,
                set[tuple[int, ...]]
            ],
            None,
            None,
        ]:
        '''Enumerate the inputs to this example as a tuple of (choice_bins, target_word)'''
        for (ignore_multiplicities, unique_bins), solution in self.solutions.items():
            yield self.choice_bins, self.target_word, ignore_multiplicities, unique_bins, solution