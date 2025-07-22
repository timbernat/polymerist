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

# convenience methods
from typing import Union
from ...genutils.trees.treeviz import treestr, Node, AbstractStyle, ContStyle 
# DEVNOTE: modules_installed() depends on anytree, so there's no need to worry about checking for it this far in this script

def available_force_fields_summary(newline : str='\n', style : Union[str, AbstractStyle]=ContStyle()) -> str:
    '''Human-readable list of all the currently-installed OpenFF ForceFields'''
    ff_dir_tree_strs : list[Node] = []
    for ff_dir_name, ff_paths in FF_PATH_REGISTRY.items():
        ff_dir_node = Node(ff_dir_name)
        for ff_path in ff_paths:
            _ = Node(ff_path.stem, parent=ff_dir_node)
        ff_dir_tree_strs.append(treestr(ff_dir_node, style=style))
        
    return newline.join(ff_dir_tree_strs)
        