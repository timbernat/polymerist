'''For assigning, transferring, and removing properties of RDKit objects'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from .atomprops import (
    atom_ids_with_prop,
    aggregate_atom_prop,
    annotate_atom_prop,
    clear_atom_props,
    RDATOM_MAGIC_PROPS,
)
from .bijection import (
    difference_rdmol,
    bijective_atom_id_iter,
)
from .rdprops import (
    isrdobj,
    detailed_rdobj_info,
    copy_rd_props,
    assign_props_from_dict,
    RDMOL_MAGIC_PROPS,
)