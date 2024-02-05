'''Tools for manipulating and extending OpenFF objects, and for interfacing with other tools and formats'''

from typing import Any
from pathlib import Path

# Locate path where OpenFF forcefields are installed
import openforcefields
FFDIR = Path(openforcefields.get_forcefield_dirs_paths()[0])


# References to all registered Toolkits
from openff.toolkit import ToolkitRegistry
from openff.toolkit import GLOBAL_TOOLKIT_REGISTRY as GTR

# Register OpenFf-compatible GNNs
## Espaloma Toolkit Wrappers
from espaloma_charge.openff_wrapper import EspalomaChargeToolkitWrapper
GTR.register_toolkit(EspalomaChargeToolkitWrapper)

## NAGL Toolkit Wrappers
from openff.nagl.toolkits import NAGLRDKitToolkitWrapper, NAGLOpenEyeToolkitWrapper
GTR.register_toolkit(NAGLRDKitToolkitWrapper)
GTR.register_toolkit(NAGLOpenEyeToolkitWrapper)

## NAGL GNN Model
from openff import nagl_models
from openff.nagl import GNNModel

NAGL_MODEL_PATH = nagl_models.list_available_nagl_models()[1] # Path(/home/timber/miniconda3/envs/polymerist-env/lib/python3.11/site-packages/openff/nagl_models/models/openff-gnn-am1bcc-0.1.0-rc.1.pt)
NAGL_MODEL_PATH = nagl_models.validate_nagl_model_path(NAGL_MODEL_PATH) # double check that this model path is still one of the valid entry point
NAGL_MODEL = GNNModel.load(NAGL_MODEL_PATH)


# Separate all registered TK Wrappers and Registries by name
TKWRAPPERS = { 
    tk_wrap.toolkit_name : tk_wrap
        for tk_wrap in GTR.registered_toolkits
}

TKREGS = {} # individually register toolkit wrappers for cases where a registry must be passed
for tk_name, tk_wrap in TKWRAPPERS.items():
    tk_reg = ToolkitRegistry()
    tk_reg.register_toolkit(tk_wrap)
    TKREGS[tk_name] = tk_reg