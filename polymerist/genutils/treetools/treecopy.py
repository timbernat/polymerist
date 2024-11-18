'''Tools for copying parts and wholes of trees, at various levels of resolution'''

from typing import Any, Optional

from functools import reduce
from anytree.node import Node
from anytree.exporter import DictExporter

from ..filters import Filter, ALWAYS_TRUE_FILTER, ALWAYS_FALSE_FILTER


def get_node_attrs(node : Node, attr_filter : Optional[Filter[str]]=None, include_name : bool=False) -> dict[str, Any]:
    '''
    Return a dict of all attributes set on a Node
    
    Parameters
    ----------
    node : Node
        An anytree.node.Node object
    attr_filter : Filter[str] (optional), default lambda x : True
        An optional criterion to decide whether an attribute should be kept
        Should be a function which accepts a single string arg and returns a bool
        Return value of True will include an attribute in the return, while value or False will exclude it
        If None, will default to always True (i.e. no attributes will be screened out)
    include_name : bool, deafult False
        Whether to include the required "name" attribute of a Node in the returned dict
        Useful to exclude when copying nodes to avoid redundancy
        By default False (i.e. "name" will be excluded from the returned attributes)

    Returns
    -------
    node_attrs : dict[str, Any]
        A dictionary keyed by attribute name whose values are the values set for repective Node attributes
    '''
    if attr_filter is None:
        attr_filter = ALWAYS_TRUE_FILTER

    if not include_name: # NOTE: need to rebind to different name to avoid recursion error
        _attr_filter = lambda attr_name : attr_filter(attr_name) and attr_name != 'name'
    else:
        _attr_filter = attr_filter

    node_attrs = {
        attr_name : attr
            for attr_name, attr in DictExporter._iter_attr_values(node) # NOTE: this is necessary to omit mangled NodeMixin info about parents and children
                if _attr_filter(attr_name)
    }
    # assert(node_copy.children == tuple())
    return node_attrs

def copy_node_unbound(node : Node, attr_filter : Optional[Filter[str]]=None) -> Node:
    '''
    Create a copy of a Node with matching attributes but no ancestors or descendents

    Parameters
    ----------
    node : Node
        An anytree.node.Node object
    attr_filter : Filter[str] (optional), default lambda x : True
        An optional criterion to decide whether an attribute should be kept
        Should be a function which accepts a single string arg and returns a bool
        Return value of True will include an attribute in the return, while value or False will exclude it
        If None, will default to always True (i.e. no attributes will be screened out)

    Returns
    -------
    node_copy : Node
        A new Node object which has the same attributes as
        the original node but none of the ancestors or descendants
    '''
    return Node(name=node.name, **get_node_attrs(node=node, attr_filter=attr_filter, include_name=False))

# NOTE: explicitly exclude a filter criterion here, as filtering (rather than stopping) may result in deleted nodes IN THE MIDDLE of a tree
def copy_tree(root : Node, stop : Optional[Filter[Node]]=None, attr_filter : Optional[Filter[str]]=None) -> Node: 
    '''
    Create a copy of an anytree Node hierarchy. Can provide filters and stop criteria to exclude nodes or whole branches
    
        Parameters
    ----------
    root : Node
        An anytree.node.Node object which is the root of a tree-like hierarchy
    stop : Filter[Node] (optional), default None
        An optional criterion to decide when to stop traversing the tree to be copied
        This criterion is not inclusive, i.e. a Node matching this criterion WILL be 
        included in the copied tree, but any child Nodes of this flagged node will not

        Should be a function which accepts a single Node arg and returns a bool
        Return value of True will exclude all subsequent nodes on a branch, while value of True will proceed with iteration and copying 
        If None, will default to always False (i.e. no extra stop conditions, full tree will be copied)
    attr_filter : Filter[str] (optional), default lambda x : True
        An optional criterion to decide whether an attribute should be kept
        Should be a function which accepts a single string arg and returns a bool
        Return value of True will include an attribute in the return, while value or False will exclude it
        If None, will default to always True (i.e. no attributes will be screened out)

    Returns
    -------
    root_new : Node
        The root Node of the copied tree structure
    '''
    if stop is None:
        stop = ALWAYS_FALSE_FILTER

    root_new = copy_node_unbound(root, attr_filter=attr_filter)
    for child in root.children: # NOTE: this also works for leaf nodes, as their "children" attrs is just an empty tuple
        if stop(child):
            continue
        child_copy = copy_tree(child, stop=stop, attr_filter=attr_filter) # recursively copy children until stop criterion
        child_copy.parent = root_new

    return root_new