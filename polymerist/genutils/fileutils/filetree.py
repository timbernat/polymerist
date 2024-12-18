'''Tools for manipulating files and directories in the file system'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Iterable

from pathlib import Path
from subprocess import Popen

from ..trees.treebase import NodeCorrespondence, compile_tree_factory
from ..trees.treeviz import treestr
from ..decorators.functional import allow_string_paths


# FILE TREES
class PathToNodeCorrespondence(NodeCorrespondence, FROMTYPE=Path):
    '''Concrete implementation of pathlib Paths as nodes in a tree'''
    def name(self, path : Path) -> str:
        return path.name
    
    def has_children(self, path : Path) -> bool:
        return path.is_dir()
    
    def children(self, path) -> Iterable[Path]:
        return path.iterdir()
    
path_tree = file_tree = allow_string_paths(
    compile_tree_factory(
        PathToNodeCorrespondence(),
        obj_attr_name='ppath' # NOTE: can't call this "path", as that clashes with an attribute of Node also called "Path"
    )
)
dir_tree = allow_string_paths(
    compile_tree_factory(
        PathToNodeCorrespondence(),
        class_alias='directory',
        obj_attr_name='dir',
        exclude_mixin=lambda path : path.is_file()
    )
)

# MODIFICATIONS TO DIRECTORIES ON DISC
def startfile(path : Path) -> None:
    '''Replacement for os.startfile() functionality, since none natively exists in Linux'''
    Popen(['xdg-open', path])

def clear_dir(path : Path) -> None:
    '''Recursively clear contents of a directory at the given path (depth-first)'''
    assert(path.is_dir())

    for sub_path in path.iterdir():
        if sub_path.is_dir():
            clear_dir(sub_path)
            sub_path.rmdir() # raises OSError if inside of target subfolder while being deleted
        else:
            sub_path.unlink()
