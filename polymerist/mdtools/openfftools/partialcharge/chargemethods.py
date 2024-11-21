'''Registry module for keeping track of which partial charging toolkit registries and related methods are available'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Type
from collections import defaultdict

from openff.toolkit.utils.base_wrapper import ToolkitWrapper
from openff.toolkit.utils.builtin_wrapper import BuiltInToolkitWrapper
from openff.toolkit.utils.rdkit_wrapper import RDKitToolkitWrapper
from openff.toolkit.utils.openeye_wrapper import OpenEyeToolkitWrapper
from openff.toolkit.utils.ambertools_wrapper import AmberToolsToolkitWrapper

from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
from espaloma_charge.openff_wrapper import EspalomaChargeToolkitWrapper

from .. import REGISTERED_TKWRAPPER_TYPES


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
