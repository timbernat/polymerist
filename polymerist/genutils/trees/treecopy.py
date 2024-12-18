'''Tools for copying parts and wholes of trees, at various levels of resolution'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Any, Optional

from anytree.node import Node
from anytree.iterators import PreOrderIter
from anytree.exporter import DictExporter

from networkx import DiGraph

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
        Should be a function which accepts a single Node arg and returns a bool
        
        This criterion is inclusive, i.e. a Node matching this criterion will NOT be 
        included in the copied tree, nor will any of its children or their children, recursively

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

def tree_to_networkx(root : Node, stop : Optional[Filter[Node]]=None, attr_filter : Optional[Filter[str]]=None) -> Node: 
    '''
    Produces a networkx.DiGraph representation on an anytree Tree

    Parameters
    ----------
    root : Node
        An anytree.node.Node object which is the root of a tree-like hierarchy
    stop : Filter[Node] (optional), default None
        An optional criterion to decide when to stop traversing the tree to be copied
        Should be a function which accepts a single Node arg and returns a bool
        
        Return value of True will exclude all subsequent nodes on a branch, while value of True will proceed with iteration and copying 
        If None, will default to always False (i.e. no extra stop conditions, full tree will be copied)

        This criterion is inclusive, i.e. a Node matching this criterion will NOT be 
        included in the copied tree, nor will any of its children or their children, recursively
    attr_filter : Filter[str] (optional), default lambda x : True
        An optional criterion to decide whether an attribute should be kept
        Should be a function which accepts a single string arg and returns a bool

        Return value of True will include an attribute in the return, while value or False will exclude it
        If None, will default to always True (i.e. no attributes will be screened out)

    Returns
    -------
    nx_tree : diGraph
        A networkx directed graph object 
    '''
    if stop is None:
        stop = ALWAYS_FALSE_FILTER

    pruned : dict[Node, bool] = {}  # which nodes should be excluded (needed for recursive branch pruning)
    node_ids : dict[Node, int] = {} # unique integer index for nodes (needed to resolve duplicate-name Nodes)
    nx_tree = DiGraph()
    for i, node in enumerate(PreOrderIter(root)): # NOTE: DON'T CHANGE ITERATOR!, implementation here relies on the topological sorting property of pre-order traversal
        node_ids[node] = i # mark node for parent identification later
        if stop(node):
            pruned[node] = True
            continue # skip over node addition once marked

        parent_idx : Optional[int] = None
        if node.parent is not None:
            parent_idx = node_ids[node.parent]
            if pruned[node.parent]: # if the parent node was pruned...
                pruned[node] = True # ...mark the current node as also pruned...
                continue            # ...and skip over its addition
        
        nx_tree.add_node(
            i, 
            name=node.name, # label=node.name # NOTE: attr "name" conflicts with pydot Node "name" init attribute; consider using "label" or similar for uniqueness?
            **get_node_attrs(node, attr_filter=attr_filter, include_name=False)
        )
        if parent_idx is not None:
            nx_tree.add_edge(parent_idx, i) # parent guaranteed to have been visited first and be mapped by pre-order topological sorting
        pruned[node] = False # only once a node is added to the graph should it be explicitly marked as not pruned

    return nx_tree