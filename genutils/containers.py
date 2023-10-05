'''Custom data containers with useful properties'''

from collections import defaultdict, Counter
from typing import Any, Iterable, TypeVar
T = TypeVar('T') # generic type variable


def RecursiveDict() -> defaultdict:
    '''Returns a defaultdict which can be recursively nested indefinitely'''
    return defaultdict(lambda : RecursiveDict())

class UnorderedRegistry:
    '''For storing and comparing unordered collections of items'''
    def __init__(self, *defaults : Iterable[T]) -> None:
        self._defaults = list(defaults) # cache defaults for future use
        self.elements : list[Counter] = []
        self.reset()

    def __repr__(self) -> str:
        elem_str = ', '.join(str(elem) for elem in self.elements)
        return f'{self.__class__.__name__}({elem_str})'

    def input_as_counter(funct): # NOTE : internal decorator for namespace, deliberately omits a "self" arg 
        '''For transmuting the first input to a method into a counter'''
        def wrapper(self, elem : Any, *args, **kwargs):
            return funct(self, Counter(elem), *args, **kwargs)
        return wrapper
    
    @input_as_counter
    def __contains__(self, elem : T) -> bool:
        '''Check for membership, compatible with Python's "in" syntax'''
        return elem in self.elements
    
    @input_as_counter
    def insert(self, elem : T) -> None:
        '''Add a new element to the registry, avoiding duplication'''
        if elem not in self.elements:
            self.elements.append(elem)
    
    @input_as_counter
    def pop(self, elem : T) -> T:
        '''Removes a registered element from the UnoderedRegistry and returns it; 
        Raises KeyError if the element was never registered in the first place'''
        idx = self.elements.index(elem)
        return self.elements.pop(idx)
    
    def clear(self) -> None:
        '''Empty all registered elements'''
        self.elements.clear()

    def reset(self) -> None:
        '''Clear contents and restore initialized defaults'''
        self.clear()
        for elem in self._defaults:
            self.insert(Counter(elem)) # TOSELF : can't use sets here, as Counter objects are mutable and therefore unhashable