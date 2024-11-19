'''Functionality for dynamically importing and inspecting Python modules and packages'''

from .pkgiter import module_hierarchy, iter_submodules, module_tree, module_tree_direct
from .pkginspect import is_package, is_module