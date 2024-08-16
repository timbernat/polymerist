'''Typehinting and generic implementations of filter (indicator) functions'''

from typing import Callable, TypeVar

T = TypeVar('T')
Filter = Callable[[T], bool] # TODO: move this to somewhere in typetools

NULL_FILTER : Filter[T] = lambda inp : False # a filter which doesn't do anything (but has the right call signature)