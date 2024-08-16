'''Generic functionality for tree-like data structures. Based on the anytree module (https://github.com/c0fec0de/anytree)'''

from typing import Any, Callable, Generic, Iterable, Optional, TypeAlias, TypeVar
from abc import ABC, abstractmethod

Filter : TypeAlias = Callable[[Any], bool] # TODO: move this to somewhere in typetools
T = TypeVar('T')

from anytree.node import Node
from anytree.exporter import DictExporter

from .decorators.classmod import register_abstract_class_attrs


def example_tree_for_tests() -> Node: # TODO: move to separate tests module eventually
    '''Produce a simplified tree for performing tests'''
    root = Node('f')
    b = Node('b', foo='bb', parent=root)
    a = Node('a', foo='aa', parent=b)
    d = Node('d', foo='dd', parent=b)
    c = Node('c', foo='cc', parent=d)
    e = Node('e', foo='ee', parent=d)
    g = Node('g', foo='gg', parent=root)
    i = Node('i', foo='ii', parent=g)
    h = Node('h', foo='hh', parent=i)

    return root


# TREE COPYING
def copy_node_isolated(node : Node) -> Node:
    '''Create a copy of a Node with only its attributes and no set ancestors or descendents'''
    node_attrs = {
        attr_name : attr
            for attr_name, attr in DictExporter._iter_attr_values(node) # NOTE: this is necessary to omit mangled NodeMixin info about parents and children
    }
    return Node(**node_attrs)

# NOTE: explicitly exclude a filter criterion here, as filtering (rather than stopping) may result in deleted nodes IN THE MIDDLE of a tree
def copy_tree(node : Node, stop : Optional[Filter]=None) -> Node: 
    '''Create a copy of an anytree Node hierarchy. Can provide filters and stop criteria to exclude nodes or whole branches'''
    if stop is None:
        stop = lambda node : False

    node_copy = copy_node_isolated(node) # make a read-only copy of JUST the current node's attributes
    assert(node_copy.children == tuple())
    if node.children is not None:
        for child in node.children:
            if stop(child):
                continue
            child_copy = copy_tree(child, stop=stop) # recursively copy children until stop criterion
            child_copy.parent = node_copy
    
    return node_copy


# INTERFACES FOR BUILDING TREES FROM OTHER CLASSES
@register_abstract_class_attrs('FROMTYPE')
class AbstractNodeCorrespondence(ABC, Generic[T]): # in concrete implementations, the type of NODETYPE should match T
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

def compile_tree_factory(node_corresp : AbstractNodeCorrespondence[T], class_alias : Optional[str]=None, obj_attr_name : Optional[str]=None) -> Callable[[T, Optional[int]], Node]:
    '''Factory method for producing a tree-generating function for the given Type''' # TODO: include blacklist
    if class_alias is None: # an alternative name to use when describing the tree creation for this class
        class_alias = node_corresp.FROMTYPE.__name__

    if obj_attr_name is None: # the name given to the Node attribute which store an instance of the given arbitrary type
        obj_attr_name = class_alias

    def compile_tree(obj : T,  max_depth : Optional[int]=None, _curr_depth : int=0) -> Node:
        # NOTE: deliberately omitting docstring here, as it will be built procedurally after defining this function
        node = Node(name=node_corresp.name(obj))
        setattr(node, obj_attr_name, obj) # keep an instance of the object directly for reference

        if node_corresp.has_children(obj) and ( # recursively add subnodes IFF
                (max_depth is None)             # 1) no depth limit is set, or
                or (_curr_depth < max_depth)    # 2) a limit IS set, but hasn't been reached yet
            ): 
            for child_obj in node_corresp.children(obj):
                sub_node = compile_tree(child_obj, max_depth=max_depth, _curr_depth=_curr_depth+1)
                sub_node.parent = node

        return node

    # annoyingly, docstrings must be string literals (this CAN'T be done inside the function definition)
    compile_tree.__doc__ = f''' 
    Compile a {class_alias} tree from a(n) {node_corresp.FROMTYPE.__name__} object

    Any sub-{class_alias} encountered will be expanded into its own tree,
    up to the specified maximum depth, or until exhaustion if max_depth=None
    '''

    return compile_tree