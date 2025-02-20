'''For assigning, transferring, and removing properties of RDKit objects'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from ._rdprops import (
    isrdobj,
    detailed_rdobj_info,
    copy_rdobj_props,
    assign_props_from_dict,
    RDMOL_MAGIC_PROPS,
    RDATOM_MAGIC_PROPS,
    RDPROP_GETTERS,
    RDPROP_SETTERS
)
