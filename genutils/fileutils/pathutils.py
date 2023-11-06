'''Utilities for editing, augmenting, and querying Paths'''

from typing import Union
from pathlib import Path
    

# PATH PROPERTY FUNCTIONS (DON'T MODIFY ANYTHING)
def _dotless(extension : str) -> str:
    '''Separate the dot from a SINGLE extension file suffix. Returns the original suffix if not dot is present'''
    return extension.split('.')[-1] 

def dotless(path : Path) -> str:
    '''Separate the dot from file path. Returns the original suffix if not dot is present'''
    return _dotless(path.suffix)

def is_empty(path : Path) -> bool:
    '''Check if a directory is empty'''
    assert(path.is_dir())
    return list(path.iterdir()) == [] # can't use "len" for generators : TODO : make this more efficient (i.e. iteration-based) for large directories


# PATH CREATION FUNCTIONS
def assemble_path(directory : Path, prefix : str, extension : str, postfix : str='') -> Path:
    '''Combine output, naming, descriptive, and filetype info to generate a complete Path'''
    if extension[0] == '.':
        extension = extension[1:] # remove leading dots if included
    path_name = f'{prefix}{"_" if postfix else ""}{postfix}.{extension}'

    return directory / path_name


# PATH CONVERSION FUNCTIONS (FOR CHANGING BETWEEN TYPES)
def aspath(path : Union[str, Path]) -> Path:
	'''Allow functions which expect Paths to also accept strings'''
	if not isinstance(path, Path):
		path = Path(path)
	return path

def asstrpath(strpath : Union[str, Path]) -> str:
	'''Allow functions which expect strings paths to also accept Paths'''
	if not isinstance(strpath, str):
		strpath = str(strpath)
	return strpath


# PATH MODIFICATION FUNCTIONS (CHANGING THE STRUCTURE OF A PATH OBJECT)
def default_suffix(path : Path, suffix : str) -> Path:
    '''Asserts that a path has a suffix, appending a specified default suffix if none exists'''
    if not path.suffix:
        path = path.with_name(f'{path.stem}.{suffix}') # ensure charge params path has correct extension

    return path

def prepend_parent(path : Path, new_parent : Path) -> Path:
    '''Prepends a parent tree to an existing path'''
    return new_parent / path

def detach_parent(path : Path, old_parent : Path) -> Path:
    '''Cuts off a parent tree from an existing path'''
    return path.relative_to(old_parent)

def exchange_parent(path : Path, old_parent : Path, new_parent : Path) -> Path:
    '''Exchanges the parent tree of a path for another parent tree'''
    return prepend_parent(path=detach_parent(path, old_parent), new_parent=new_parent)

def local_rename(path : Path, new_name : str) -> Path:
    '''Performs file rename relative to the parent directory (NOT the cwd)'''
    return path.rename(path.with_name(new_name))

def local_restem(path : Path, new_stem : str) -> Path:
    '''Performs file rename relative to the parent directory (NOT the cwd), preserving the extension of the original file'''
    return path.rename(path.with_stem(new_stem))
