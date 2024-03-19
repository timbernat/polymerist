'''Tools for manipulating and extending OpenFF objects, and for interfacing with other tools and formats'''

from typing import Any
from pathlib import Path

import openforcefields
from openff.toolkit import ToolkitRegistry
from openff.toolkit import GLOBAL_TOOLKIT_REGISTRY as GTR
from openff.toolkit.utils.base_wrapper import ToolkitWrapper
from openff.toolkit.utils.utils import all_subclasses
from openff.toolkit.utils.exceptions import LicenseError, ToolkitUnavailableException

from openff.toolkit.utils.openeye_wrapper import OpenEyeToolkitWrapper
from espaloma_charge.openff_wrapper import EspalomaChargeToolkitWrapper
from openff.nagl.toolkits import NAGLRDKitToolkitWrapper, NAGLOpenEyeToolkitWrapper


# FORCE FIELD AND ToolkitWrapper REFERENCE
FFDIR = Path(openforcefields.get_forcefield_dirs_paths()[0]) # Locate path where OpenFF forcefields are installed
ALL_IMPORTABLE_TKWRAPPERS = all_subclasses(ToolkitWrapper) # References to every registered ToolkitWrapper and ToolkitRegistry 

# CHECKING FOR OpenEye
try:
    _ = OpenEyeToolkitWrapper()
    _OE_TKWRAPPER_IS_AVAILABLE = True
    OEUnavailableException = None
except (LicenseError, ToolkitUnavailableException) as error:
    _OE_TKWRAPPER_IS_AVAILABLE = False
    OEUnavailableException = error # catch and record relevant error message for use (rather than trying to replicate it elsewhere)

# Register OpenFF-compatible GNN ToolkitWrappers
GTR.register_toolkit(EspalomaChargeToolkitWrapper)
GTR.register_toolkit(NAGLRDKitToolkitWrapper)
if _OE_TKWRAPPER_IS_AVAILABLE:
    GTR.register_toolkit(NAGLOpenEyeToolkitWrapper)


# GENERATE LOOKUP DICTS FOR EVERY REGISTERED ToolkitWrappers and ToolkitRegistry
REGISTERED_TKWRAPPER_TYPES = [type(tkwrapper) for tkwrapper in GTR.registered_toolkits]
TKWRAPPERS = { # NOTE : this must be done AFTER any new registrations to thr GTR (e.g. after registering GNN ToolkitWrappers)
    tk_wrap.toolkit_name : tk_wrap
        for tk_wrap in GTR.registered_toolkits
}
TKREGS = {} # individually register toolkit wrappers for cases where a registry must be passed
for tk_name, tk_wrap in TKWRAPPERS.items():
    tk_reg = ToolkitRegistry()
    tk_reg.register_toolkit(tk_wrap)
    TKREGS[tk_name] = tk_reg