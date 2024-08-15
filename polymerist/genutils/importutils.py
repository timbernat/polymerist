'''Utilities for automating submodule import and logger creation'''

import logging
LOGGER = logging.getLogger(__name__)

from types import ModuleType
from typing import Generator, Iterable, Optional

import pkgutil
import importlib

from enum import StrEnum
from itertools import chain


# FILETREE CHARACTERS
_TREE_WIDTH : int = 4
assert(_TREE_WIDTH > 0)

class TreeChars(StrEnum):
    '''Box characters for representing connections between components of a hierarchcal structure'''
    SPACE   = ' '*_TREE_WIDTH
    DASH    = '\u2500'*_TREE_WIDTH
    PIPE    = '\u2502' + ' '*(_TREE_WIDTH - 1)
    BRANCH  = '\u251C' + '\u2500'*(_TREE_WIDTH - 1)
    ELBOW   = '\u2514' + '\u2500'*(_TREE_WIDTH - 1)
    ELBOW_R = '\u2514' + '\u2500'*(_TREE_WIDTH - 1) # direction-agnostic alias
    ELBOW_L = '\u2500'*(_TREE_WIDTH - 1) + '\u2518'


# BASE SUBMODULE GENERATORS
def module_by_pkg_str(pkg_str : str) -> ModuleType:
    '''Load module from dot-separated module name string (intended to be called on __package__ attr)'''
    return importlib.import_module(pkg_str)

def iter_submodule_info(module : ModuleType, recursive : bool=True, blacklist : Optional[Iterable[str]]=None) -> Generator[tuple[ModuleType, str, bool], None, None]:
    '''Generate all submodules of a given module, yielding a tuple of (the module, the module name, and whether the module is also a package).
    If the "recursive" flag is set, will generate ALL possible submodules in the tree recursively'''
    if blacklist is None:
        blacklist = []
    
    for _loader, submodule_name, submodule_ispkg in pkgutil.iter_modules(module.__path__):
        if submodule_name in blacklist:
            continue # skip over import blacklisted modules

        try:
            submodule = importlib.import_module(f'{module.__package__}.{submodule_name}')
        except ModuleNotFoundError:
            continue

        yield (submodule, submodule_name, submodule_ispkg)
        if submodule_ispkg and recursive:
            yield from iter_submodule_info(submodule, recursive=True, blacklist=blacklist)

def iter_submodules(module : ModuleType, recursive : bool=True, blacklist : Optional[Iterable[str]]=None) -> Generator[ModuleType, None, None]:
    '''Generate all submodules of a given module. If the "recursive" flag is set, will generate ALL possible submodules in the tree recursively'''
    for submodule, *_submodule_info in iter_submodule_info(module, recursive=recursive, blacklist=blacklist):
        yield submodule # only yield submodule; for more compact iteration


# TOOLS FOR REGISTERING AND PULLING INFO FROM SUBMODULES
def register_submodules(module : ModuleType, recursive : bool=True, blacklist : Optional[Iterable[str]]=None) -> None:
    '''Registers submodules of a given module into it's own namespace (i.e. autoimports submodules)'''
    for (submodule, submodule_name, submodule_ispkg) in iter_submodule_info(module, recursive=False, blacklist=blacklist): # initially only iterate on one level to keep track of parent module
        setattr(module, submodule_name, submodule)
        if submodule_ispkg and recursive:
            register_submodules(submodule, recursive=recursive, blacklist=blacklist)

def _module_hierarchy(module : ModuleType, recursive : bool=True, blacklist : Optional[Iterable[str]]=None, _prefix : str='') -> Generator[str, None, None]:
    '''Returns an iterable of level strings representing a module hierarchy, with each level of nesting indicated by a series of pipes and dashes'''
    end_sentinel = (object(), object(), object()) # used to unambiguously check whether the initial generator is in fact empty
    module_iter = chain(iter_submodule_info(module, recursive=False, blacklist=blacklist), [end_sentinel])

    module_info = next(module_iter) # primer for hierarchy iteration, need to first check 
    reached_end = bool(module_info == end_sentinel) # if module hierarchy is empty, don't bother trying to iterate
    while not reached_end:
        (submodule, submodule_name, submodule_ispkg) = module_info # unpacking to keep current module_info values in namespace (and for convenience)
        module_info = next(module_iter) # peek ahead to check if the current module is in fact the last at this layer
        if module_info == end_sentinel:
            splitter, extension = TreeChars.ELBOW, TreeChars.SPACE
            reached_end = True
        else:
            splitter, extension = TreeChars.BRANCH, TreeChars.PIPE
        
        yield _prefix + splitter + submodule_name
        if submodule_ispkg and recursive:
            yield from _module_hierarchy(submodule, recursive=recursive, blacklist=blacklist, _prefix=_prefix + extension) # retrieve partial output

def module_hierarchy(module : ModuleType, recursive : bool=True, blacklist : Optional[Iterable[str]]=None) -> str:
    return '\n'.join(_module_hierarchy(module, recursive=recursive, blacklist=blacklist))

def submodule_loggers(module : ModuleType, recursive : bool=True, blacklist : Optional[Iterable[str]]=None) -> dict[str, Optional[logging.Logger]]:
    '''Produce a dict of any Logger objects present in each submodule. Can optionally generate recursively and blacklist certain modules'''
    return {
        submodule_name : getattr(submodule, 'LOGGER', None) # default to None rather than raising Exception
            for (submodule, submodule_name, submodule_ispkg) in iter_submodule_info(module, recursive=recursive, blacklist=blacklist) 
    }