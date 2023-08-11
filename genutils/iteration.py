'''Tools for simplifying iteration over collections of items'''

from typing import Any, Callable, Generator, Iterable, TypeVar, Union

from itertools import islice, product as cartesian_product
from collections import deque

from .decorators.functional import optional_in_place


# PROPERTIES OF ITERABLES
def iter_len(itera : Iterable):
    '''Get size of an iterable object where ordinary len() call is invalid (namely a generator)
    Note that this will "use up" a generator upon iteration'''
    return sum(1 for _  in itera)


# READ-ONLY ITERATION
T = TypeVar('T') # generic type of object(s) in iterator
T1 = TypeVar('T1') # generic type of object(s) in first iterator
T2 = TypeVar('T2') # generic type of object(s) in second iterator

def sliding_window(items : Iterable[T], n : int=1) -> Generator[tuple[T], None, None]:
    '''Generates sliding windows of width n over an iterable collection of items
    E.g. : sliding_window('ABCDE', 3) --> (A, B, C), (B, C, D), (C, D, E)
    '''
    it = iter(items)
    window = deque(islice(it, n), maxlen=n)
    if len(window) == n:
        yield tuple(window)
    for x in it:
        window.append(x)
        yield tuple(window)

def swappable_loop_order(iter1 : Iterable[T1], iter2 : Iterable[T2], swap : bool=False) -> Union[Iterable[tuple[T1, T2]], Iterable[tuple[T2, T1]]]:
    '''Enables dynamic swapping of the order of execution of a 2-nested for loop'''
    order = [iter1, iter2] if not swap else [iter2, iter1]
    for pair in cartesian_product(*order):
        yield pair[::(-1)**swap] # reverse order of pair (preserves argument identity)

def progress_iter(itera : Iterable[T], key : Callable[[T], str]=lambda x : x) -> Iterable[tuple[str, T]]:
    '''Iterate through'''
    N = len(itera) # TODO : extend this to work for generators / consumables
    for i, item in enumerate(itera):
        yield (f'{key(item)} ({i + 1} / {N})', item) # +1 converts to more human-readable 1-index for step count


# READ-WRITE ITERATION (MODIFIES THE OBJECTS BEING ITERATED OVER)
K = TypeVar('K') # generic type for dict key
V = TypeVar('V') # generic type for dict value
O = TypeVar('O') # generic type for an object passed to a function

def asiterable(arg_val : Union[O, Iterable[O]]) -> Iterable[O]:
	'''Permits functions expecting iterable arguments to accept singular values'''
	if not isinstance(arg_val, Iterable):
		arg_val = (arg_val,) # turn into single-item tuple (better for memory)
	return arg_val

@optional_in_place
def modify_dict(path_dict : dict[K, V], modifier_fn : Callable[[K, V], Any]) -> None:
    '''Recursively modifies all values in a dict in-place according to some function'''
    for key, val in path_dict.items():
        if isinstance(val, dict): # recursive call if sub-values are also dicts with Paths
            modify_dict(val, modifier_fn)
        else:
            path_dict[key] = modifier_fn(key, val) 

def sort_dict_by_values(targ_dict : dict, reverse : bool=False) -> dict[Any, Any]:
    '''Sort a dictionary according to the values of each key'''
    return { # sort dict in ascending order by size
        key : targ_dict[key]
            for key in sorted(targ_dict, key=lambda k : targ_dict[k], reverse=reverse)
    }