'''Utilities for generating, labelling, editing, and transforming RDKit molecules and other RDKit objects'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from .rdkdraw import (
    set_rdkdraw_size,
    enable_substruct_highlights,
    disable_substruct_highlights,
    enable_kekulized_drawing,
    disable_kekulized_drawing,
)

# Define default drawing settings - TODO: find way to configure these at package level on import/install
set_rdkdraw_size(300, aspect=3/2)
enable_substruct_highlights()
disable_kekulized_drawing() # explicitly show when rings are aromatic or Kekule on drawing