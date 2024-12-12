'''Extensions, interfaces, and convenience methods built around the functionality in the OpenFF software stack'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

# Subpackage-wide precheck to see if OpenFF is even usable in the first place
from ...genutils.importutils.dependencies import modules_installed, MissingPrerequisitePackage
if not modules_installed('openff', 'openff.toolkit'):
    raise MissingPrerequisitePackage(
        importing_package_name=__spec__.name,
        use_case='OpenFF addons',
        install_link='https://docs.openforcefield.org/projects/toolkit/en/stable/installation.html',
        dependency_name='openff-toolkit',
        dependency_name_formal='the OpenFF software stack',
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