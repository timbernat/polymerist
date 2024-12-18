'''For representing, generating, and modifying information about groups of monomer fragments'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from .repr import MonomerGroup # make monomer representation available at the module level
from .fragments import (
    PE_FRAGMENTS,
    MPD_TMC_FRAGMENTS,
    PEG_PLGA_FRAGMENTS,
)