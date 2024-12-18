'''Generic functionality for tree-like data structures. Based on the anytree module (https://github.com/c0fec0de/anytree)'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from .treebase import NodeCorrespondence, compile_tree_factory
from .treecopy import get_node_attrs, copy_tree, tree_to_networkx
from .treeviz import treestr