'''Unit tests for tree-related functionality'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Any
import pytest

from anytree.node import Node
from anytree.search import find as find_node
from anytree.iterators import AbstractIter, PreOrderIter, PostOrderIter, LevelOrderIter

from networkx import dfs_preorder_nodes

from polymerist.genutils.trees.treecopy import get_node_attrs, copy_node_unbound, copy_tree, tree_to_networkx


@pytest.fixture
def example_tree() -> Node: 
    '''Produce a simplified tree for performing tests'''
    root = Node('f')
    b = Node('b', parent=root, foo='bb')
    a = Node('a', parent=b,    foo='aa', bar='spam')
    d = Node('d', parent=b,    foo='dd')
    c = Node('c', parent=d,    foo='cc', bar='spam')
    e = Node('e', parent=d,    foo='ee', baz=[1,2,3])
    g = Node('g', parent=root, foo='gg', bar='spam', baz=[1,2,3])
    i = Node('i', parent=g,    foo='ii')
    h = Node('h', parent=i,    foo='hh')
    a_dup = Node('a', foo='a+a', parent=g) # testing how nodes with duplicate names are handled

    return root

@pytest.fixture
def g_node(example_tree) -> Node:
    '''Returns a particularly attribute-rich and neighbor-rich node for testing'''
    return find_node(example_tree, filter_=lambda node : node.name == 'g')

@pytest.fixture
def node_name() -> str:
    '''Allows for configuration of test name values (test should pass with arbitrary names)'''
    return 'test'

@pytest.fixture
def node_attrs() -> dict[str, Any]:
    '''Allows for configuration of test attr values (test should pass with arbitrary attrs)'''
    return {
        'foo'  : 'bar',
        'baz'  : 42,
        'spam' : 'eggs',
    }

def test_get_node_attrs_no_name(node_name : str, node_attrs : dict[str, Any]) -> None:
    '''Test that correct Node attributes are fetched (excluding Node name)'''
    node = Node(name=node_name, **node_attrs)
    attrs = get_node_attrs(node, include_name=False)
    assert attrs == node_attrs

def test_get_node_attrs_with_name(node_name : str, node_attrs : dict[str, Any]) -> None:
    '''Test that correct Node attributes are fetched (including Node name)'''
    node = Node(name=node_name, **node_attrs)
    attrs = get_node_attrs(node, include_name=True)
    assert attrs == {'name' : node_name, **node_attrs}

def test_get_node_attrs_attr_name_filter(node_name : str, node_attrs : dict[str, Any]) -> None:
    '''Test that filter conditions are correctly applied when fetching Node attributes'''
    node = Node(name=node_name, **node_attrs)
    attrs = get_node_attrs(node, attr_filter=lambda attr_name : attr_name != 'foo', include_name=False)
    assert 'foo' not in attrs


def test_copy_node_unbound_values(g_node : Node) -> None:
    '''Test that copy_node_unbound() correctly copies Node attributes'''
    copied_node = copy_node_unbound(g_node)
    assert get_node_attrs(copied_node) == get_node_attrs(g_node)

def test_copy_node_unbound_relatives(g_node : Node) -> None:
    '''Test that copy_node_unbound() correctly removes ancestors and children from node'''
    copied_node = copy_node_unbound(g_node)
    assert copied_node.parent is None and not copied_node.children


@pytest.mark.parametrize('iter_type', [PreOrderIter, PostOrderIter, LevelOrderIter])
def test_copy_tree(example_tree : Node, iter_type : AbstractIter) -> None:
    '''Test that trees structures are exactly copied with no filters'''
    copied_tree = copy_tree(example_tree, stop=None, attr_filter=None)
    assert all(
        get_node_attrs(node_orig) == get_node_attrs(node_copied) 
            for node_orig, node_copied in zip(iter_type(example_tree), iter_type(copied_tree))
    )

@pytest.mark.parametrize('iter_type', [PreOrderIter, PostOrderIter, LevelOrderIter])
def test_copy_tree_mod_attr(example_tree : Node, iter_type : AbstractIter) -> None:
    '''Test that modifying attributes on copied trees does NOT affect original'''
    copied_tree = copy_tree(example_tree, stop=None, attr_filter=None)
    TARG_ATTR : str = 'foo'
    for node in iter_type(copied_tree): # modify attributes
        if hasattr(node, TARG_ATTR):
            setattr(node, TARG_ATTR, getattr(node, TARG_ATTR) + '__')

    assert all( # attributes should NOT be equal, since they should've only been changed on the copy
        getattr(node_orig, TARG_ATTR) != getattr(node_copied, TARG_ATTR)
            for node_orig, node_copied in zip(iter_type(example_tree), iter_type(copied_tree))
                if hasattr(node_orig, TARG_ATTR)
    )

def test_copy_tree_stop(example_tree : Node) -> None:
    '''Test that stop conditions for tree copying are respected and targetted branches are pruned'''
    copied_tree = copy_tree(example_tree, stop=lambda node : node.name == 'a', attr_filter=None)
    assert all(
        copied_node.name != 'a'
            for copied_node in PreOrderIter(copied_tree) # in this case, the iteration order uniquely doesn't matter, only care that all nodes are traversed 
    )


def test_tree_to_networkx(example_tree : Node) -> None:
    '''Test that conversion to networkx.DiGraph faithfully reproduces node order and attributes'''
    nxtree = tree_to_networkx(example_tree)
    assert all(
        nxtree.nodes[i] == get_node_attrs(at_node, include_name=True)
            for i, at_node in zip(dfs_preorder_nodes(nxtree), PreOrderIter(example_tree)) # preorder should match, regardless of implementation
    )