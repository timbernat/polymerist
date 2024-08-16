'''Tools for copying parts and wholes of trees, at various levels of resolution'''

from typing import Optional

from anytree.node import Node
from anytree.exporter import DictExporter

from ..filters import Filter, NULL_FILTER


def copy_node_attrs(node : Node, attr_filter : Filter[str]=NULL_FILTER) -> Node:
    '''Create a copy of a Node with only its attributes and no set ancestors or descendents'''
    node_attrs = {
        attr_name : attr
            for attr_name, attr in DictExporter._iter_attr_values(node) # NOTE: this is necessary to omit mangled NodeMixin info about parents and children
                if attr_filter(attr_name)
    }
    # assert(node_copy.children == tuple())
    return Node(**node_attrs)

# NOTE: explicitly exclude a filter criterion here, as filtering (rather than stopping) may result in deleted nodes IN THE MIDDLE of a tree
def copy_tree(node : Node, stop : Optional[Filter[Node]]=NULL_FILTER, attr_filter : Filter[str]=NULL_FILTER) -> Node: 
    '''Create a copy of an anytree Node hierarchy. Can provide filters and stop criteria to exclude nodes or whole branches'''
    node_copy = copy_node_attrs(node, attr_filter=attr_filter) # make a read-only copy of JUST the current node's attributes
    for child in node.children: # NOTE: this also works for leaf nodes, as their "children" attrs is just an empty tuple
        if stop(child):
            continue
        child_copy = copy_tree(child, stop=stop) # recursively copy children until stop criterion
        child_copy.parent = node_copy

    return node_copy