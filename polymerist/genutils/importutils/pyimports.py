'''For inspecting and managing toplevel imports within Python files and modules'''

from typing import Optional
from types import ModuleType
from dataclasses import dataclass, field

import ast
from pathlib import Path

from .decorators.functional import allow_string_paths # TODO: see if there's anyway to bypass a relative import here


@dataclass
class ImportedObjectInfo:
    '''For encapsulating info about an object imported in a Python file'''
    object_name   : str
    object_alias  : Optional[str]  = field(default=None)
    parent_module : Optional[str]  = field(default=None)
    source_file   : Optional[Path] = field(default=None)
    line_number   : Optional[int]  = field(default=None)
    is_relative   : Optional[bool] = field(default=None)

@allow_string_paths
def extract_imports_from_pyfile(pyfile_path : Path) -> list[ImportedObjectInfo]:
    '''Compiles info from all Python imports in a Python (.py) file'''
    if not pyfile_path.is_file():
        raise ValueError('Cannot interpret non-file path as file')
    
    with pyfile_path.open('r') as pyfile:
        module_root = ast.parse(pyfile.read(), pyfile_path)

    import_info : list[ImportedObjectInfo] = []
    for syntax_node in module_root.body:
        if not isinstance(syntax_node, (ast.Import, ast.ImportFrom)):
            continue
        
        if isinstance(syntax_node, ast.ImportFrom):
            parent_module : Optional[str] = syntax_node.module
            is_relative   : bool = (syntax_node.level != 0)
        else: # implcictly, this MUST be an ast.Import instance by the initial exclusion check
            parent_module : Optional[str] = None
            is_relative   : bool = False

        for imported_object in syntax_node.names:
            import_info.append(
                ImportedObjectInfo(
                    object_name=imported_object.name,
                    object_alias=imported_object.asname,
                    parent_module=parent_module,
                    source_file=pyfile_path,
                    line_number=syntax_node.lineno, # .end_lineno
                    is_relative=is_relative,
                )
            )
    
    return import_info

@allow_string_paths
def extract_imports_from_dir(source_dir : Path) -> list[ImportedObjectInfo]:
    '''Compiles info from all Python imports in any and all Python files in a directory'''
    if not source_dir.is_dir():
        raise ValueError('Cannot interpret non-directory path as directory')
    
    import_infos = []
    for pyfile_path in source_dir.glob('**/*.py'):
        import_infos.extend(extract_imports_from_pyfile(pyfile_path))

    return import_infos

def extract_imports_from_module(module : ModuleType) -> list[ImportedObjectInfo]:
    '''Compiles info from all Python imports in a Python (.py) file'''
    # TODO: find more reliable/canonical way to tell packages and bare modules/scripts apart
    if hasattr(module, '__file__') and getattr(module, '__file__') is not None:
        return extract_imports_from_pyfile(module.__file__)

    if hasattr(module, '__path__') and getattr(module, '__path__') is not None:
        return extract_imports_from_dir(module.__path__[0])