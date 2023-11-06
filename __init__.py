'''Polymer-Oriented Library Yielding Structure Assignment, Calculation of CHARges, Interchange, and Data Elucidation (second revision)'''

__author__ = 'Timotej Bernat'
__version__ = '2.0'

from .genutils.importutils import register_submodules, module_by_pkg_str

_MODULE_SELF = module_by_pkg_str(__package__) # keep reference to own module
register_submodules(_MODULE_SELF, recursive=True, blacklist=['analysis'])