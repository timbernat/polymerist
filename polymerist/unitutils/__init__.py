'''Physical constants, dimensional analysis, and unit conversion utilities'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from .dimensions import(
    MissingUnitsError,
    hasunits,
    strip_units,
    is_volume,
    quantities_approx_equal,
)
