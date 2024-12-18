'''Typehinting and generic implementations of filter (indicator) functions'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Callable, TypeVar

T = TypeVar('T')
Filter = Callable[[T], bool] # TODO: move this to somewhere in typetools

# "TRIVIAL" filters to use as defaults when the Filter call signature is required
ALWAYS_TRUE_FILTER  : Filter[T] = lambda inp : True 
ALWAYS_FALSE_FILTER : Filter[T] = lambda inp : False

# aliases for trivial filters for many use-cases
NEVER_FALSE_FILTER = MARK_ALL_FILTER  = FLAG_ALL_FILTER =  ALWAYS_TRUE_FILTER
NEVER_TRUE_FILTER  = MARK_NONE_FILTER = FLAG_NONE_FILTER = ALWAYS_FALSE_FILTER