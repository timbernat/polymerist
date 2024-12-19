'''Tools for iterating over and extracting information from Python package hierarchies'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import logging
LOGGER = logging.getLogger(__name__)

from types import ModuleType
from typing import Generator, Iterable, Optional, Union

from importlib import import_module
from pkgutil import iter_modules as _iter_modules

from anytree.node import Node
from anytree.render import AbstractStyle, ContStyle
from anytree.iterators import PreOrderIter

from .pkginspect import module_stem, is_package
from .dependencies import MissingPrerequisitePackage
from ..trees.treebase import NodeCorrespondence, compile_tree_factory
from ..trees.treeviz import treestr



# HIERARCHICAL MODULE TREE GENERATION
class ModuleToNodeCorrespondence(NodeCorrespondence, FROMTYPE=ModuleType):
    '''Concrete implementation of Python modules and packages as nodes in a tree'''
    def name(self, module : ModuleType) -> str:
        return module_stem(module) # TODO: find Python package-spec compliant way of extracting this easily
    
    def has_children(self, module : ModuleType) -> bool:
        return is_package(module)
    
    def children(self, module : ModuleType) -> Iterable[ModuleType]:
        for _loader, module_name, ispkg in _iter_modules(module.__path__, prefix=module.__name__+'.'):
            try:
                submodule = import_module(module_name)
                yield submodule
            except (ModuleNotFoundError, MissingPrerequisitePackage):
                continue

module_tree = compile_tree_factory(
    ModuleToNodeCorrespondence(),
    class_alias='package',
    obj_attr_name='module',
    exclude_mixin=lambda module : module_stem(module).startswith('_'),
)

# BACKWARDS-COMPATIBLE PORTS OF LEGACY IMPORTUTILS FUNCTIONS
def module_tree_direct(module : ModuleType, recursive : bool=True, blacklist : Optional[Iterable[str]]=None) -> Node:
    '''Produce a tree from the Python package hierarchy starting with a given module
    
    Parameters
    ----------
    module : ModuleType
        The "root" module to begin importing from
        Represented in the Node object returned by this function
    recursive : bool, default=True
        Whether or not to recursively import modules from subpackages and add them to the tree
    blacklist : list[str] (optional), default None
        List of module names to exclude from tree building
        If provided, will exclude any modules whose names occur in this list

    Returns
    -------
    modtree : Node
        The root node of the module tree, corresponding to the module object passed to "module"
    '''
    if blacklist is None:
        blacklist = []

    return module_tree(
        module,
        max_depth=None if recursive else 1,
        exclude=lambda module : module_stem(module) in blacklist,
    )

def iter_submodules(module : ModuleType, recursive : bool=True, blacklist : Optional[Iterable[str]]=None) -> Generator[ModuleType, None, None]:
    '''
    Generates all modules which can be imported from the given toplevel module
    
    Parameters
    ----------
    module : ModuleType
        The "root" module to begin importing from
        Represented in the Node object returned by this function
    recursive : bool, default=True
        Whether or not to recursively import modules from subpackages and add them to the tree
    blacklist : list[str] (optional), default None
        List of module names to exclude from tree building
        If provided, will exclude any modules whose names occur in this list

    Returns
    -------
    submodules : Generator[ModuleType]
        A generator which yields modules in traversal pre-order as they appear wihin the package hierarchy
    '''
    modtree = module_tree_direct(module, recursive=recursive, blacklist=blacklist)
    for module_node in PreOrderIter(modtree):
        yield module_node.module

def iter_submodule_info(module : ModuleType, recursive : bool=True, blacklist : Optional[Iterable[str]]=None) -> Generator[tuple[ModuleType, str, bool], None, None]:
    '''
    Generates information about all modules which can be imported from the given toplevel module
    Namely, yields the module object, module name, and whether or not the module is a package

    Parameters
    ----------
    module : ModuleType
        The "root" module to begin importing from
        Represented in the Node object returned by this function
    recursive : bool, default=True
        Whether or not to recursively import modules from subpackages and add them to the tree
    blacklist : list[str] (optional), default None
        List of module names to exclude from tree building
        If provided, will exclude any modules whose names occur in this list

    Returns
    -------
    submodule_info : Generator[ModuleType, str, bool]
        A generator which yields modules info in traversal pre-order as they appear wihin the package hierarchy
        yields 3-tuples containing ModuleType objects, module names, and whether the current module is also a subpackage
    '''
    modtree = module_tree_direct(module, recursive=recursive, blacklist=blacklist)
    for module_node in PreOrderIter(modtree):
        yield module_node.module, module_node.name, module_node.is_leaf

def register_submodules(module : ModuleType, recursive : bool=True, blacklist : Optional[Iterable[str]]=None) -> None:
    '''
    Registers all submodules of a given module into it's own namespace (i.e. autoimports submodules)
    
    Parameters
    ----------
    module : ModuleType
        The "root" module to begin importing from
        Represented in the Node object returned by this function
    recursive : bool, default=True
        Whether or not to recursively import modules from subpackages and add them to the tree
    blacklist : list[str] (optional), default None
        List of module names to exclude from tree building
        If provided, will exclude any modules whose names occur in this list

    Returns
    -------
    None
    '''
    for submodule in iter_submodules(module, recursive=recursive, blacklist=blacklist):
        setattr(module, submodule.__name__, submodule)

def module_hierarchy(module : ModuleType, recursive : bool=True, blacklist : Optional[Iterable[str]]=None, style : Union[str, AbstractStyle]=ContStyle()) -> str:
    '''
    Generates a printable string which summarizes a Python packages hierarchy. Reminiscent of GNU tree output

    Parameters
    ----------
    module : ModuleType
        The "root" module to begin importing from
        Represented in the Node object returned by this function
    recursive : bool, default=True
        Whether or not to recursively import modules from subpackages and add them to the tree
    blacklist : list[str] (optional), default None
        List of module names to exclude from tree building
        If provided, will exclude any modules whose names occur in this list
    style : str or AbstractStyle
        An element drawing style for the final tree structure printout

    Returns
    -------
    module_summary : str
        Printable string which displays the package structure
    '''
    modtree = module_tree_direct(module, recursive=recursive, blacklist=blacklist)
    return treestr(modtree, style=style)

