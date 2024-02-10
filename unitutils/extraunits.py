'''Defining units which, for one reason or another, are not defined in Pint or OpenMM units'''

import sys
_MODULE = sys.modules[__name__]

from openmm.unit import Unit, ScaledUnit
from openmm.unit.prefix import define_prefixed_units


# ENERGY
from openmm.unit import joule
from scipy.constants import electron_volt as electron_volt_in_joules

electronvolt_base_unit = ScaledUnit(electron_volt_in_joules, joule, 'electronvolt', 'eV')
electronvolt = eV = Unit({electronvolt_base_unit : 1.0})
define_prefixed_units(electronvolt_base_unit, module=_MODULE)


# VISCOSITY
from openmm.unit import centimeter, gram, second

poise_base_unit = ScaledUnit(1.0, gram/(centimeter*second), 'poise', 'P')
poise = Unit({poise_base_unit : 1.0})
define_prefixed_units(poise_base_unit, module=_MODULE)