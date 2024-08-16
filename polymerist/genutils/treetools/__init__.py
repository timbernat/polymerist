'''Generic functionality for tree-like data structures. Based on the anytree module (https://github.com/c0fec0de/anytree)'''

from typing import Any, Callable, TypeVar

T = TypeVar('T')
Filter = Callable[[T], bool] # TODO: move this to somewhere in typetools
NULL_FILTER : Filter[T] = lambda inp : False # a filter which doesn't do anything (but has the right call signature)
