'''Utilities for automating submodule import and logger creation'''

from types import ModuleType
from typing import Generator, Iterable, Optional

import pkgutil
import importlib

import logging


# BASE SUBMODULE GENERATORS
def module_by_pkg_str(pkg_str : str) -> ModuleType:
    '''Load module from dot-separated module name string intended to be called on __package__ attr)'''
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
    for submodule, _submodule_info in iter_submodule_info(module, recursive=recursive, blacklist=blacklist):
        yield submodule # only yield submodule; for more compact iteration


# TOOLS FOR REGISTERING AND PUUL INFO FROM SUBMODULES
def register_submodules(module : ModuleType, recursive : bool=True, blacklist : Optional[Iterable[str]]=None) -> None:
    '''Registers submodules of a given module into it's own namespace (i.e. autoimports submodules)'''
    for (submodule, submodule_name, submodule_ispkg) in iter_submodule_info(module, recursive=False, blacklist=blacklist): # initially only iterate on one level to keep track of parent module
        setattr(module, submodule_name, submodule)
        if submodule_ispkg and recursive:
            register_submodules(submodule, recursive=recursive, blacklist=blacklist)

def module_hierarchy(module : ModuleType, recursive : bool=True, blacklist : Optional[Iterable[str]]=None, _indent_level : int=0, _string_agg : Optional[list[str]]=None) -> str:
    '''Returns a string representing a module hierarchy, with each level of indentation representing one more level of nesting'''
    if _string_agg is None:
        _string_agg = []
    tabs = _indent_level*'\t'

    for (submodule, submodule_name, submodule_ispkg) in iter_submodule_info(module, recursive=False, blacklist=blacklist): # initially only iterate on one level to keep track of parent module
        _string_agg.append(f'{tabs}{submodule_name}')
        if submodule_ispkg and recursive:
            _partial = module_hierarchy(submodule, recursive=recursive, blacklist=blacklist, _indent_level=_indent_level + 1, _string_agg=_string_agg) # retrieve partial output

    return '\n'.join(_string_agg)

def submodule_loggers(module : ModuleType, recursive : bool=True, blacklist : Optional[Iterable[str]]=None) -> dict[str, Optional[logging.Logger]]:
    '''Produce a dict of any Logger objects present in each submodule. Can optionally generate recursively and blacklist certain modules'''
    return {
        submodule_name : getattr(submodule, 'LOGGER', None) # default to None rather than raising Exception
            for (submodule, submodule_name, submodule_ispkg) in iter_submodule_info(module, recursive=recursive, blacklist=blacklist) 
    }