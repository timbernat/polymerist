'''References to global OpenFF objects and files'''

from typing import Any
from pathlib import Path

# Locate path where OpenFF forcefields are installed
import openforcefields
FFDIR = Path(openforcefields.get_forcefield_dirs_paths()[0])


# References to all registered Toolkits
from openff.toolkit import ToolkitRegistry
from openff.toolkit import GLOBAL_TOOLKIT_REGISTRY as GTR

## Register GNN Toolkits
from espaloma_charge.openff_wrapper import EspalomaChargeToolkitWrapper
GTR.register_toolkit(EspalomaChargeToolkitWrapper)

# from openff.nagl.toolkits import NAGLRDKitToolkitWrapper, NAGLOpenEyeToolkitWrapper
# GTR.register_toolkit(NAGLRDKitToolkitWrapper)
# GTR.register_toolkit(NAGLOpenEyeToolkitWrapper)

## Extract name-keyed TK wrappers and Registries
TKWRAPPERS = { 
    tk_wrap.toolkit_name : tk_wrap
        for tk_wrap in GTR.registered_toolkits
}

TKREGS = {} # individually register toolkit wrappers for cases where a registry must be passed
for tk_name, tk_wrap in TKWRAPPERS.items():
    tk_reg = ToolkitRegistry()
    tk_reg.register_toolkit(tk_wrap)
    TKREGS[tk_name] = tk_reg

# def __getattr__(name : str) -> Any:
#     '''Hook for defining property-like attributes of this module'''
#     if name == 'TOOLKITS':
#         return { # cast as lambda to return updated version in case toolskits change
#             tk.toolkit_name : tk
#                 for tk in GTR.registered_toolkits
#         }
    
#     raise AttributeError(f"Module '{__name__}' has no attribute '{name}'")
