'''Unit tests for trees'''

from anytree.node import Node

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