'''General-purpose utilities related to SMILES and SMARTS string manipulations'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from .cleanup import (
    # validation
    Smiles,
    is_valid_SMILES,
    Smarts,
    is_valid_SMARTS,
    # custom Exceptions
    InvalidChemicalLineNotation,
    InvalidSMILES,
    InvalidSMARTS,
    InvalidInChI,
    # canonicalization
    canonical_SMILES_from_mol,
    expanded_SMILES,
)