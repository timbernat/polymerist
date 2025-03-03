'''Protocols for interfaces not explicitly covered by typing or collections.abc in the stdlib'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import (
    # Generics
    Any,
    Generic,
    TypeVar,
    # Protocols
    Iterable,
    Protocol,
    runtime_checkable,
)
ReturnType = TypeVar('ReturnType')


@runtime_checkable
class IndexableIterable(Protocol, Generic[ReturnType]):
    '''Any container-like class which supports indexing and iteration'''
    def __getitem__(key : Any) -> ReturnType:
        ...
        
    def __iter__(self) -> Iterable[ReturnType]:
        ...