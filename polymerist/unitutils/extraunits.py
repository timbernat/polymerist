'''Defining units which, for one reason or another, are not defined in Pint or OpenMM units'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import sys
_MODULE = sys.modules[__name__]

from math import pi
from openmm.unit import Unit, BaseUnit, ScaledUnit
from openmm.unit.prefix import define_prefixed_units

from openmm.unit import second_base_unit, single_item_amount_base_unit
from openmm.unit import coulomb_base_unit, ampere_base_unit, watt_base_unit
from openmm.unit import time_dimension, charge_dimension, amount_dimension
from openmm.unit import joule, erg, hartree
from openmm.unit import centimeter, gram, second
from openmm.unit import meter, candela

from scipy.constants import (
    c as speed_of_light,
    hbar as reduced_planck_const,
    electron_volt as electron_volt_in_joules,
)


# PREFIX EXPANSIONS NOT FOUND BY DEFAULT IN OPENMM UNITS
define_prefixed_units(coulomb_base_unit, module=_MODULE)
define_prefixed_units(ampere_base_unit , module=_MODULE)
define_prefixed_units(watt_base_unit   , module=_MODULE)


# UNITS FOR AMOUNTS OF ITEMS
percent_base_unit = BaseUnit(amount_dimension, 'percent', '%')
percent_base_unit.define_conversion_factor_to(single_item_amount_base_unit, 1 / 100)
percent = Unit({percent_base_unit : 1})

dozen_base_unit = BaseUnit(amount_dimension, 'dozen', 'doz')
dozen_base_unit.define_conversion_factor_to(single_item_amount_base_unit, 12)
dozen = Unit({dozen_base_unit : 1})


# PHYSICIST units - TODO : make atomic_unit_system UnitSystem subclass
## ENERGY
electronvolt_base_unit = ScaledUnit(electron_volt_in_joules, joule, 'electronvolt', 'eV')
electronvolt = eV = Unit({electronvolt_base_unit : 1.0})
define_prefixed_units(electronvolt_base_unit, module=_MODULE)

## TIME
atomic_time_base_unit = BaseUnit(time_dimension, 'atomic_time', 'au')
atomic_time_base_unit.define_conversion_factor_to(second_base_unit, reduced_planck_const / hartree.conversion_factor_to(joule))
atomic_time = Unit({atomic_time_base_unit : 1.0})

## FREQUENCY
hertz_base_unit = ScaledUnit(1.0, second**-1, 'hertz', 'Hz')
hertz = Unit({hertz_base_unit : 1.0})
define_prefixed_units(hertz_base_unit, module=_MODULE) # TODO : add conversion to radian per second (base units are incompatible though!!)


# PHOTOMETRIC UNITS
lumen_base_unit = ScaledUnit(1/(4*pi), candela, 'lumen', 'lm')
lumen = lumens = Unit({lumen_base_unit : 1})
# deliberately omitted prefix registration for now

lux_base_unit = ScaledUnit(1.0, lumen / meter**2, 'lux', 'lx')
lux = Unit({lux_base_unit : 1.0})
# deliberately omitted prefix registration for now


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