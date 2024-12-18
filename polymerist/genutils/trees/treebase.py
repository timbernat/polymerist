'''Interfaces for encoding arbitrary classes into tree-like data structures'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Any, Callable, Generic, Iterable, Optional, TypeAlias, TypeVar
from abc import ABC, abstractmethod

from anytree.node import Node

from ..decorators.classmod import register_abstract_class_attrs
from ..filters import Filter, ALWAYS_FALSE_FILTER

T = TypeVar('T')


@register_abstract_class_attrs('FROMTYPE') # TODO: figure out way to parameterize Generic T here with the type passed as FROMTYPE
class NodeCorrespondence(ABC, Generic[T]): 
    '''Abstract base for implementing how to build an anytree Node tree for an arbitrary class'''
    @abstractmethod
    def name(self, obj : T) -> str:
        '''Define how to obtain a string name'''
        pass

    @abstractmethod
    def has_children(self, obj : T) -> bool:
        '''Define how to check if an object can produce children in the first place before attempting to do so'''
        pass

    @abstractmethod
    def children(self, obj : T) -> Optional[Iterable[T]]:
        '''Define how to obtain node children from an instance
        Should return NoneType if the instance is "leaf-like"'''
        pass

def compile_tree_factory(
        node_corresp : NodeCorrespondence[T],
        class_alias : Optional[str]=None,
        obj_attr_name : Optional[str]=None,
        exclude_mixin : Optional[Filter[T]]=None,
    ) -> Callable[[T, Optional[int], Optional[Filter[T]]], Node]:
    '''
    Factory method for producing a tree-generating function from a NodeCorrespondence
    
    Parameters
    ----------
    node_corresp : NodeCorrespondence[T]
        Definition of a correpondence between an arbitrary type and a Tree Node
    class_alias : str (optional)
        Name of the corresponding class to inject into docstring
        If not provided, will default to the __name__ of the class wrapped by node_corresp
    obj_attr_name : str (optional)
        The name of the Node attribute to which a copy of the class instance should be bound
        If not provided, will default to the value of class_alias
    exclude_mixin : Filter[T] (optional)
        An optional "master" filter to mix into any 

    Returns
    -------
    compile_tree : Callable[[T, Optional[int], Optional[Filter[T]]], Node]
        Factory function which takes an instance of type T, builds a Tree from it, and returns the root Node
    '''
    node_type_name = node_corresp.FROMTYPE.__name__
    if class_alias is None: # an alternative name to use when describing the tree creation for this class
        class_alias = node_type_name

    if obj_attr_name is None: # the name given to the Node attribute which store an instance of the given arbitrary type
        obj_attr_name = class_alias
    if hasattr(Node, obj_attr_name):
        raise AttributeError(f'Invalid value for obj_attr_name; attribute "{obj_attr_name}" clashes with existing attribute Node.{obj_attr_name}')

    if exclude_mixin is None:
        exclude_mixin = ALWAYS_FALSE_FILTER

    def compile_tree(
        obj : node_corresp.FROMTYPE,
        max_depth : Optional[int]=None,
        exclude : Filter[node_corresp.FROMTYPE]=ALWAYS_FALSE_FILTER,
        _curr_depth : int=0
    ) -> Node:
        # NOTE: deliberately omitting docstring here, as it will be built procedurally after defining this function
        node = Node(name=node_corresp.name(obj))
        setattr(node, obj_attr_name, obj) # keep an instance of the object directly for reference

        if node_corresp.has_children(obj) and ( # recursively add subnodes IFF
                (max_depth is None)             # 1) no depth limit is set, or
                or (_curr_depth < max_depth)    # 2) a limit IS set, but hasn't been reached yet
            ): 
            for child_obj in node_corresp.children(obj):
                if not (exclude(child_obj) or exclude_mixin(child_obj)):
                    sub_node = compile_tree(child_obj, max_depth=max_depth, exclude=exclude, _curr_depth=_curr_depth+1)
                    sub_node.parent = node

        return node

    # annoyingly, docstrings must be string literals (this CAN'T be done inside the function definition)
    compile_tree.__doc__ = f''' 
    Compile a {class_alias} tree from a(n) {node_type_name} object

    Any sub-{class_alias} encountered will be expanded into its own tree,
    up to the specified maximum depth, or until exhaustion if max_depth=None

    Parameters
    ----------
    obj : {node_type_name}
        A(n) instance of a {class_alias}
    max_depth : int (optional)
        Maximum allowed height of a constructed tree from the root
        If None (as default), no limit is set
    exclude : Filter[{node_type_name}] (optional)
        An optional filter function to reduce the size of the constructed tree
        Must accept a(n) {node_type_name} instance as single argument
        Should return True when a node is to be excluded, and False otherwise
    '''

    return compile_tree