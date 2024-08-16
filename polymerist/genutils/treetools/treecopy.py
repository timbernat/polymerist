'''Tools for copying parts and wholes of trees, at various levels of resolution'''

from typing import Optional

from anytree.node import Node
from anytree.exporter import DictExporter

from . import Filter, NULL_FILTER


def copy_node_isolated(node : Node) -> Node:
    '''Create a copy of a Node with only its attributes and no set ancestors or descendents'''
    node_attrs = {
        attr_name : attr
            for attr_name, attr in DictExporter._iter_attr_values(node) # NOTE: this is necessary to omit mangled NodeMixin info about parents and children
    }
    return Node(**node_attrs)

# NOTE: explicitly exclude a filter criterion here, as filtering (rather than stopping) may result in deleted nodes IN THE MIDDLE of a tree
def copy_tree(node : Node, stop : Optional[Filter[Node]]=None) -> Node: 
    '''Create a copy of an anytree Node hierarchy. Can provide filters and stop criteria to exclude nodes or whole branches'''
    if stop is None:
        stop = NULL_FILTER

    node_copy = copy_node_isolated(node) # make a read-only copy of JUST the current node's attributes
    assert(node_copy.children == tuple())
    if node.children is not None:
        for child in node.children:
            if stop(child):
                continue
            child_copy = copy_tree(child, stop=stop) # recursively copy children until stop criterion
            child_copy.parent = node_copy
    
    return node_copy