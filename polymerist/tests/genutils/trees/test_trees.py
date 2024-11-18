'''Unit tests for trees'''

from typing import Any
import pytest

from anytree.node import Node
from anytree.search import find as find_node

from polymerist.genutils.trees.treecopy import get_node_attrs, copy_node_unbound, copy_tree


@pytest.fixture
def example_tree() -> Node: 
    '''Produce a simplified tree for performing tests'''
    root = Node('f')
    b = Node('b', parent=root, foo='bb')
    a = Node('a', parent=b,    foo='aa', bar='bar')
    d = Node('d', parent=b,    foo='dd')
    c = Node('c', parent=d,    foo='cc', bar='bar')
    e = Node('e', parent=d,    foo='ee', baz=[1,2,3])
    g = Node('g', parent=root, foo='gg', bar='bar', baz=[1,2,3])
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
    copy_node = copy_node_unbound(g_node)
    assert get_node_attrs(copy_node) == get_node_attrs(g_node)

def test_copy_node_unbound_relatives(g_node : Node) -> None:
    '''Test that copy_node_unbound() correctly removes ancestors and children from node'''
    copy_node = copy_node_unbound(g_node)
    assert copy_node.parent is None and not copy_node.children