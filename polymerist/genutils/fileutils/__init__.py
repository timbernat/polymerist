'''Utilities for manipulating Path-like objects and interfacing with directories and files'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from .pathutils import (
    is_empty,
    assemble_path,
    allow_string_paths,
    allow_pathlib_paths,
    prepend_parent,
    detach_parent,
    exchange_parent,
    local_rename,
    local_restem,
)
from .filetree import (
    file_tree,
    dir_tree,
    startfile,
    temporary_cd,
    clear_dir,
)