'''Tools for interfacing and representing arbitrary external classes with tree-like data structures'''

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
    ) -> Callable[[T, Optional[int]], Node]:
    '''Factory method for producing a tree-generating function for the given Type''' # TODO: include blacklist
    if class_alias is None: # an alternative name to use when describing the tree creation for this class
        class_alias = node_corresp.FROMTYPE.__name__

    if obj_attr_name is None: # the name given to the Node attribute which store an instance of the given arbitrary type
        obj_attr_name = class_alias

    def compile_tree(obj : node_corresp.FROMTYPE,  max_depth : Optional[int]=None, exclude : Filter[node_corresp.FROMTYPE]=ALWAYS_FALSE_FILTER, _curr_depth : int=0) -> Node:
        # NOTE: deliberately omitting docstring here, as it will be built procedurally after defining this function
        node = Node(name=node_corresp.name(obj))
        setattr(node, obj_attr_name, obj) # keep an instance of the object directly for reference

        if node_corresp.has_children(obj) and ( # recursively add subnodes IFF
                (max_depth is None)             # 1) no depth limit is set, or
                or (_curr_depth < max_depth)    # 2) a limit IS set, but hasn't been reached yet
            ): 
            for child_obj in node_corresp.children(obj):
                if not exclude(child_obj):
                    sub_node = compile_tree(child_obj, max_depth=max_depth, exclude=exclude, _curr_depth=_curr_depth+1)
                    sub_node.parent = node

        return node

    # annoyingly, docstrings must be string literals (this CAN'T be done inside the function definition)
    compile_tree.__doc__ = f''' 
    Compile a {class_alias} tree from a(n) {node_corresp.FROMTYPE.__name__} object

    Any sub-{class_alias} encountered will be expanded into its own tree,
    up to the specified maximum depth, or until exhaustion if max_depth=None
    '''

    return compile_tree