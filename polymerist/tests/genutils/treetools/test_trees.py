'''Unit tests for trees'''

from typing import Any
import pytest

from anytree.node import Node

from polymerist.genutils.treetools import treecopy


@pytest.fixture
def example_tree() -> Node: 
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
    a_dup = Node('a', foo='a+a', parent=g) # testing how nodes with duplicate names are handled

    return root

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
    attrs = treecopy.get_node_attrs(node, include_name=False)
    assert attrs == node_attrs

def test_get_node_attrs_with_name(node_name : str, node_attrs : dict[str, Any]) -> None:
    '''Test that correct Node attributes are fetched (including Node name)'''
    node = Node(name=node_name, **node_attrs)
    attrs = treecopy.get_node_attrs(node, include_name=True)
    assert attrs == {'name' : node_name, **node_attrs}

def test_get_node_attrs_attr_name_filter(node_name : str, node_attrs : dict[str, Any]) -> None:
    '''Test that correct Node attributes are fetched (including Node name)'''
    node = Node(name=node_name, **node_attrs)
    attrs = treecopy.get_node_attrs(node, attr_filter=lambda attr_name : attr_name != 'foo', include_name=False)
    assert 'foo' not in attrs