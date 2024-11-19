'''Generic functionality for tree-like data structures. Based on the anytree module (https://github.com/c0fec0de/anytree)'''

from .treebase import NodeCorrespondence, compile_tree_factory
from .treecopy import get_node_attrs, copy_tree, tree_to_networkx
from .treeviz import treestr