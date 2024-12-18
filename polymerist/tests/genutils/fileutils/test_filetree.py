'''Unit tests for filetree operations'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import pytest
from pathlib import Path
from anytree.iterators import PreOrderIter, LevelOrderGroupIter

from polymerist.genutils.importutils.pkginspect import get_dir_path_within_package
from polymerist.tests import data as testdata

from polymerist.genutils.fileutils.filetree import path_tree, dir_tree


DUMMY_DIR_INFO : dict[int, bool] = { # expected depth and dir status of test directory
    'dummy_dir' : (0, True),
    'subdir1'   : (1, True),
    'bar.txt'   : (1, False),
    'foo.dat'   : (1, False),
    'subdir2'   : (2, True),
    'spam.dat'  : (2, False),
    'baz.txt'   : (3, False),
}

@pytest.fixture
def dummy_dir_path() -> Path:
    return get_dir_path_within_package('dummy_dir', testdata)

@pytest.mark.parametrize('depth', [0, 1, 2, 3, 4])
def test_path_tree_depth(dummy_dir_path, depth : int) -> None:
    '''Test that max_depth restrictions are correctly applied'''
    tree = path_tree(dummy_dir_path, max_depth=depth)
    assert tree.height <= depth

@pytest.mark.parametrize('depth', [0, 1, 2, 3, 4])
def test_path_tree_output(dummy_dir_path, depth : int) -> None:
    '''Test that file names and path/dir status is correctly recorded'''
    tree = path_tree(dummy_dir_path, max_depth=depth)
    for node_group in LevelOrderGroupIter(tree):
        for node in node_group:
            expected = DUMMY_DIR_INFO[node.name]
            encountered = (node.depth, node.ppath.is_dir())
            assert(expected == encountered)

def test_path_tree_exclude(dummy_dir_path) -> None:
    '''Check that path_tree correctly excludes paths by filter conditions'''
    tree = path_tree(dummy_dir_path, max_depth=None, exclude=lambda path : path.suffix == '.dat')
    for node in PreOrderIter(tree):
        assert node.ppath.suffix != '.dat'

def test_dir_tree(dummy_dir_path) -> None:
    '''Test that dir_tree correctly filters down to just directories'''
    tree = dir_tree(dummy_dir_path, max_depth=None)
    for node in PreOrderIter(tree):
        assert node.dir.is_dir() # note the necessary change of attribute from "ppath" to "dir"

def test_dir_tree_exclude(dummy_dir_path) -> None:
    '''Test that dir_tree correctly excludes directories by filter conditions'''
    tree = dir_tree(dummy_dir_path, max_depth=None, exclude=lambda dir : dir.name.startswith('subdir'))
    for node in PreOrderIter(tree):
        assert not node.name.startswith('subdir') # note the necessary change of attribute from "ppath" to "dir"