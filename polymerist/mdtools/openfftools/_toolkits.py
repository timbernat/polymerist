'''For dynamically determining and cataloging which ToolkitWrappers (and accompanying functionality) are available'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

# Core OpenFF toolkit component registration
from typing import Union
from collections import defaultdict

from openff.toolkit.utils.utils import all_subclasses
from openff.toolkit.utils.base_wrapper import ToolkitWrapper
from openff.toolkit.utils.toolkit_registry import ToolkitRegistry 
from openff.toolkit.utils.toolkits import (
    OPENEYE_AVAILABLE,
    RDKIT_AVAILABLE,
    AMBERTOOLS_AVAILABLE,
    GLOBAL_TOOLKIT_REGISTRY,
)
GTR = GLOBAL_TOOLKIT_REGISTRY # alias for brevity

from ...genutils.importutils.dependencies import modules_installed
_REGISTER_TOOLKITS_TO_GLOBAL : bool = True # TODO: find way to avoid setting this config parameter directly in code


# Helper functions
def toolkit_wrapper_is_registered(toolkit_wrapper : Union[ToolkitWrapper, type[ToolkitWrapper]], toolkit_registry : ToolkitRegistry=GTR) -> bool:
    '''Check whether a ToolkitRegistry instance has already registered a given ToolkitWrapper subclass'''
    if not isinstance(toolkit_wrapper, type):   # ToolkitWrapper TYPES are needed for this check; any instances therefore...
        toolkit_wrapper = type(toolkit_wrapper) # ...will have their respective types extracted
    if not issubclass(toolkit_wrapper, ToolkitWrapper):
        raise TypeError(f'Expected a ToolkitWrapper instance or subclass, instead received object of type {toolkit_wrapper.__name__}')
        
    return any(isinstance(tkwrapper, toolkit_wrapper) for tkwrapper in toolkit_registry.registered_toolkits)

# Setup of initial containers for OpenFF module info
ALL_IMPORTABLE_TKWRAPPERS : list[type[ToolkitWrapper]] = all_subclasses(ToolkitWrapper) # NOTE: just because you can import these, doesn't necessarily mean they can be instantiated
ALL_AVAILABLE_TKWRAPPERS  : list[type[ToolkitWrapper]] = []
CHARGE_METHODS_BY_TOOLKIT : dict[type[ToolkitWrapper], list[str]] = defaultdict(list)

# Toolkit-specific registrations which depend on available packages
## BuiltIn (not particularly useful in and of itself, but nice to know it's accessible)
if modules_installed('openff.toolkit'): # this check is idempotent to initial OpenFF check, but is nice to have for consistency between all ToolkitWrappers below
    from openff.toolkit.utils.builtin_wrapper import BuiltInToolkitWrapper
    
    ALL_AVAILABLE_TKWRAPPERS.append(BuiltInToolkitWrapper)
    CHARGE_METHODS_BY_TOOLKIT[BuiltInToolkitWrapper] = [
        charge_method
            for charge_method in BuiltInToolkitWrapper._supported_charge_methods
    ]
## RDKit
if modules_installed('rdkit') and RDKIT_AVAILABLE:
    from openff.toolkit.utils.rdkit_wrapper import RDKitToolkitWrapper
    
    ALL_AVAILABLE_TKWRAPPERS.append(RDKitToolkitWrapper)
    CHARGE_METHODS_BY_TOOLKIT[RDKitToolkitWrapper] = [
        charge_method
            for charge_method in RDKitToolkitWrapper._supported_charge_methods
    ]
## Ambertools
if modules_installed('pdb4amber') and AMBERTOOLS_AVAILABLE: # turns out "ambertools" can't actually be imported as a module, need to check for peripheral modules which are better behaved instead
    from openff.toolkit.utils.ambertools_wrapper import AmberToolsToolkitWrapper
    
    ALL_AVAILABLE_TKWRAPPERS.append(AmberToolsToolkitWrapper)
    CHARGE_METHODS_BY_TOOLKIT[AmberToolsToolkitWrapper] = [
        charge_method
            for charge_method in AmberToolsToolkitWrapper._supported_charge_methods
    ]
## OpenEye
if modules_installed('openeye.oechem', 'openeye.oeomega') and OPENEYE_AVAILABLE:
    from openff.toolkit.utils.openeye_wrapper import OpenEyeToolkitWrapper
    
    ALL_AVAILABLE_TKWRAPPERS.append(OpenEyeToolkitWrapper)
    CHARGE_METHODS_BY_TOOLKIT[OpenEyeToolkitWrapper] = [
        charge_method
            for charge_method in OpenEyeToolkitWrapper._supported_charge_methods
    ]
## NAGL - extracting available charge methods is a little different for GNN toolkits
if modules_installed('openff.nagl', 'openff.nagl_models'):
    from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
    
    ALL_AVAILABLE_TKWRAPPERS.append(NAGLToolkitWrapper)
    CHARGE_METHODS_BY_TOOLKIT[NAGLToolkitWrapper] = [
        model_path.name # need to extract dynamically from Paths
            for model_path in NAGLToolkitWrapper.list_available_nagl_models()
    ]
## Espaloma - extracting available charge methods is a little different for GNN toolkits
if modules_installed('espaloma_charge'):
    from espaloma_charge.openff_wrapper import EspalomaChargeToolkitWrapper
    
    ALL_AVAILABLE_TKWRAPPERS.append(EspalomaChargeToolkitWrapper)
    CHARGE_METHODS_BY_TOOLKIT[EspalomaChargeToolkitWrapper] = [
        'espaloma-am1bcc'
    ] # this is, at this of writing, the only available method for EspalomaCharge and unfortunately not accessible dynamically
    
# Post-registration info compilation
## Compiling name-based lookups for available ToolkitWrappers and registering all to a local
POLYMERIST_TOOLKIT_REGISTRY = ToolkitRegistry() # retain a local registry separate from GLOBAL_TOOLKIT_REGISTRY
TKWRAPPERS      : dict[str,      ToolkitWrapper ] = {} 
TKWRAPPER_TYPES : dict[str, type[ToolkitWrapper]] = {} 

for tkwrapper_type in ALL_AVAILABLE_TKWRAPPERS:
    tkwrapper_instance = tkwrapper_type() # instantiate toolkit wrapper class
    POLYMERIST_TOOLKIT_REGISTRY.register_toolkit(tkwrapper_instance)
    # if requested, also mirror all found toolkits to the Global ToolkitRegistry
    if _REGISTER_TOOLKITS_TO_GLOBAL and not toolkit_wrapper_is_registered(tkwrapper_type, GLOBAL_TOOLKIT_REGISTRY): # make registration idempotent
        GLOBAL_TOOLKIT_REGISTRY.register_toolkit(tkwrapper_instance)
    
    # register to name-based lookup dict
    TKWRAPPERS[     tkwrapper_type._toolkit_name] = tkwrapper_instance
    TKWRAPPER_TYPES[tkwrapper_type._toolkit_name] = tkwrapper_type

## Compiling registry of which partial charge methods are supported by which toolkits
TOOLKITS_BY_CHARGE_METHOD : dict[str, list[type[ToolkitWrapper]]] = defaultdict(list) # also compile inverse mapping (compiled once available toolkits are known)
for tkwrapper_type, supported_methods in CHARGE_METHODS_BY_TOOLKIT.items():
    if (tkwrapper_type in ALL_AVAILABLE_TKWRAPPERS): # exclude non-registered toolkits to avoid confusion
        for method in supported_methods:
            TOOLKITS_BY_CHARGE_METHOD[method].append(tkwrapper_type)
TOOLKITS_BY_CHARGE_METHOD = dict(TOOLKITS_BY_CHARGE_METHOD) # convert to pure dict for typing purposes
