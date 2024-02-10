'''Defining units which, for one reason or another, are not defined in Pint or OpenMM units'''

import sys
_MODULE = sys.modules[__name__]

from openmm.unit import Unit, BaseUnit, ScaledUnit
from openmm.unit.prefix import define_prefixed_units

from openmm.unit import second_base_unit, coulomb_base_unit
from openmm.unit import time_dimension, charge_dimension
from openmm.unit import joule, erg, hartree
from openmm.unit import centimeter, gram, second
from scipy.constants import (
    c as speed_of_light,
    hbar as reduced_planck_const,
    electron_volt as electron_volt_in_joules,
)


# PREFIX EXPANSIONS NOT FOUND BY DEFAULT IN OPENMM UNITS
define_prefixed_units(coulomb_base_unit, module=_MODULE)

# PHYSICIST units - TODO : make atomic_unit_system UnitSystem subclass
## ENERGY
electronvolt_base_unit = ScaledUnit(electron_volt_in_joules, joule, 'electronvolt', 'eV')
electronvolt = eV = Unit({electronvolt_base_unit : 1.0})
define_prefixed_units(electronvolt_base_unit, module=_MODULE)

## TIME
atomic_time_base_unit = BaseUnit(time_dimension, 'atomic_time', 'au')
atomic_time_base_unit.define_conversion_factor_to(second_base_unit, reduced_planck_const / hartree.conversion_factor_to(joule))
atomic_time = Unit({atomic_time_base_unit : 1.0})


# CGS units
## VISCOSITY
poise_base_unit = ScaledUnit(1.0, gram/(centimeter*second), 'poise', 'P')
poise = Unit({poise_base_unit : 1.0})
define_prefixed_units(poise_base_unit, module=_MODULE)

## CHARGE
statcoulomb_base_unit = BaseUnit(charge_dimension, 'statcoulomb', 'statC') # TOSELF : evidently, this correspondence isn;t always exact
statcoulomb_base_unit.define_conversion_factor_to(coulomb_base_unit, 1/(10*speed_of_light)) # derived from the speed of light
statcoulomb = Unit({statcoulomb_base_unit : 1.0})
# deliberately omitted prefix registration for now

# ELECTRICAL POTENTIAL

statvolt_base_unit = ScaledUnit(1.0, erg / statcoulomb, 'statvolt', 'statV')
statvolt = Unit({statvolt_base_unit : 1.0})
# deliberately omitted prefix registration for now