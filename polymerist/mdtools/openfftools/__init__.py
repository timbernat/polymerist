'''Extensions, interfaces, and convenience methods built around the functionality in the OpenFF software stack'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

# Subpackage-wide precheck to see if OpenFF is even usable in the first place
from ...genutils.importutils.dependencies import modules_installed
if not modules_installed('openff', 'openff.toolkit'):
    raise ModuleNotFoundError(
        f'''
        OpenFF packages which are required to utilitize {__name__} not found in current environment
        Please follow installation instructions at https://docs.openforcefield.org/projects/toolkit/en/stable/installation.html, then retry import
        '''
    )
    
# Import of toplevel OpenFF object registries
from ._forcefields import (
    FFDIR,
    FF_DIR_REGISTRY,
    FF_PATH_REGISTRY,
)
from ._toolkits import (
    ## toolkit registries
    GLOBAL_TOOLKIT_REGISTRY, GTR,
    POLYMERIST_TOOLKIT_REGISTRY,
    ## catalogues of available toolkit wrappers
    ALL_IMPORTABLE_TKWRAPPERS,
    ALL_AVAILABLE_TKWRAPPERS,
    TKWRAPPERS,
    TKWRAPPER_TYPES,
    ## registry of partial charge methods by 
    CHARGE_METHODS_BY_TOOLKIT,
    TOOLKITS_BY_CHARGE_METHOD,
)