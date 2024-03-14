'''Registry module for keeping track of which partial charging toolkit registries and related methods are available'''

from typing import Type
from collections import defaultdict

from openff.toolkit.utils.base_wrapper import ToolkitWrapper
from openff.toolkit.utils.rdkit_wrapper import RDKitToolkitWrapper
from openff.toolkit.utils.builtin_wrapper import BuiltInToolkitWrapper
from openff.toolkit.utils.openeye_wrapper import OpenEyeToolkitWrapper
from openff.toolkit.utils.ambertools_wrapper import AmberToolsToolkitWrapper

from openff import nagl_models
from openff.nagl import GNNModel
from espaloma_charge.openff_wrapper import EspalomaChargeToolkitWrapper

from polymerist.openfftools import REGISTERED_TKWRAPPER_TYPES


# REFERENCE MAPPING BETWEEN PARTIAL CHARGE METHODS AND SUPPORTING TOOLKIT WRAPPERS
SUPPORTED_PARTIAL_CHARGE_METHODS_BY_TOOLKIT : dict[Type[ToolkitWrapper], list[str]]= { # TOSELF : this unfortunately cannot be accessed dynamically as a attribute of each TollkitWrapper class
    BuiltInToolkitWrapper : [
        'zeros',
        'formal_charge',
    ],
    RDKitToolkitWrapper : [
        'mmff94',
        'gasteiger'
    ],
    AmberToolsToolkitWrapper : [
        'am1bcc',
        'am1-mulliken',
        'gasteiger',
    ],
    OpenEyeToolkitWrapper : [
        'am1bcc',
        'am1-mulliken',
        'gasteiger',
        'mmff94',
        'am1bccnosymspt',
        'am1elf10',
        'am1bccelf10',
    ],
    EspalomaChargeToolkitWrapper : [
        'espaloma-am1bcc'
    ]
}

TOOLKITS_BY_CHARGE_METHOD : dict[str, list[Type[ToolkitWrapper]]] = defaultdict(list)
for tkwrapper_type, supported_methods in SUPPORTED_PARTIAL_CHARGE_METHODS_BY_TOOLKIT.items():
    if (tkwrapper_type in REGISTERED_TKWRAPPER_TYPES): # exclude non-registered toolkits to avoid confusion
        for method in supported_methods:
            TOOLKITS_BY_CHARGE_METHOD[method].append(tkwrapper_type)


## NAGL GNN Model
NAGL_MODEL_PATH = nagl_models.list_available_nagl_models()[1] # Path(/home/timber/miniconda3/envs/polymerist-env/lib/python3.11/site-packages/openff/nagl_models/models/openff-gnn-am1bcc-0.1.0-rc.1.pt)
NAGL_MODEL_PATH = nagl_models.validate_nagl_model_path(NAGL_MODEL_PATH) # double check that this model path is still one of the valid entry point
NAGL_MODEL = GNNModel.load(NAGL_MODEL_PATH)