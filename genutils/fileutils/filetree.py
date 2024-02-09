'''Tools for manipulating files and directories in the file system'''

from typing import Callable, Optional

from pathlib import Path
from subprocess import Popen


# Path filetree functions (act on file system and directories)
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
