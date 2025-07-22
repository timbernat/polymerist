'''Tools for simplifying iteration over collections of items'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Any, Callable, Generator, Hashable, Iterable, TypeVar, Union

from operator import mul
from functools import reduce
from collections import deque
from itertools import (
    islice,
    combinations,
    product as cartesian_product,
)
from .decorators.functional import optional_in_place


# READ-ONLY ITERATION
T = TypeVar('T') # generic type of object(s) in iterator
T1 = TypeVar('T1') # generic type of object(s) in first iterator
T2 = TypeVar('T2') # generic type of object(s) in second iterator

def iter_len(itera : Iterable):
    '''Get size of an iterable object where ordinary len() call is invalid (namely a generator)
    Note that this will "use up" a generator upon iteration'''
    return sum(1 for _  in itera)

def product(itera : Iterable[T]) -> T:
    '''Multiplicative analogue to builtin sum()'''
    return reduce(mul, itera) # assumes that the type T implements __mul__/__rmul__

def sliding_window(items : Iterable[T], n : int=1) -> Generator[tuple[T], None, None]:
    '''Generates sliding windows of width n over an iterable collection of items
    E.g. : sliding_window('ABCDE', 3) --> (A, B, C), (B, C, D), (C, D, E)
    '''
    it = iter(items)
    window = deque(islice(it, n), maxlen=n)
    if len(window) == n:
        yield tuple(window)

    for x in it: # implicit else
        window.append(x)
        yield tuple(window)

def powerset(items : Iterable[T], exclude_empty : bool=False, exclude_full : bool=False) -> Generator[tuple[T], None, None]:
    '''Generate all possible subsets of a set. Can optionally exclude the empty set or the complete set (or both)'''
    base_set = set(items)
    for i in range(int(exclude_empty), len(base_set) + int(not exclude_full)):
        yield from combinations(base_set, i)

def asiterable(arg_val : Union[T, Iterable[T]]) -> Iterable[T]:
	'''Permits functions expecting iterable arguments to accept singular values'''
	if not isinstance(arg_val, Iterable):
		arg_val = (arg_val,) # turn into single-item tuple (better for memory)
	return arg_val

# DICTIONARY ITERATION
V = TypeVar('V') # generic type for dict value

@optional_in_place
def modify_dict(some_dict : dict[Hashable, V], modifier_fn : Callable[[Hashable, V], Any]) -> None:
    '''Recursively modifies all values in a dict in-place according to some function'''
    for key, value in some_dict.items():
        if isinstance(value, dict): # recursive call if sub-values are also dicts with Paths
            modify_dict(value, modifier_fn)
        else:
            some_dict[key] = modifier_fn(key, value) 

def sort_dict_by_keys(targ_dict : dict, reverse : bool=False) -> dict[Any, Any]:
    '''Sort a dictionary according to the values of each key'''
    return { # sort dict in ascending order by size
        key : targ_dict[key]
            for key in sorted(targ_dict, reverse=reverse)
    }

def sort_dict_by_values(targ_dict : dict, reverse : bool=False) -> dict[Any, Any]:
    '''Sort a dictionary according to the values of each key'''
    return { # sort dict in ascending order by size
        key : targ_dict[key]
            for key in sorted(targ_dict, key=lambda k : targ_dict[k], reverse=reverse)
    }
    
def cartesian_grid(param_dict : dict[str, Iterable[Any]]) -> Generator[dict[str, Any], None, None]:
    '''
    Accepts a dict maping keys to all valid parameter values for that key
    Generates all possible choices of parameters as dicts with exactly one parameter value per key
    '''
    for param_choices in cartesian_product(*param_dict.values()):
        yield {
            param_name : param_value
                for param_name, param_value in zip(param_dict.keys(), param_choices)
        }