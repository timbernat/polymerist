'''Reference for LAMMPS unit styles, as listed in https://docs.lammps.org/units.html'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import ClassVar, Optional
from dataclasses import dataclass, field
from ...genutils.decorators.classmod import register_subclasses

# sensible units
from openmm.unit import *
from openmm.unit import Unit, Quantity
# NOTE : these refixed imports are redundant, and are just there to stop linters from whinging
from openmm.unit import kilogram, attogram, picogram                     # mass
from openmm.unit import centimeter, micrometer, nanometer                # length
from openmm.unit import microsecond, nanosecond, picosecond, femtosecond # time

from ...unitutils.extraunits import picocoulomb
from ...unitutils.extraunits import electronvolt, atomic_time    # physicist units
from ...unitutils.extraunits import statcoulomb, statvolt, poise # CGS units


# Base container class for LAMMPS unit styles
@dataclass(frozen=True) # should to be immutable to avoid mistakes
@register_subclasses(key_attr='STYLE_NAME')
class LAMMPSUnitStyle:
    '''For encapsulating which units a particular LAMMPS unit style uses'''
    STYLE_NAME : ClassVar[str] = 'BASE'

    mass              : Unit = field(init=False) # user should never be able to initialize any of these
    distance          : Unit = field(init=False)
    time              : Unit = field(init=False)
    energy            : Unit = field(init=False)
    velocity          : Unit = field(init=False)
    force             : Unit = field(init=False)
    torque            : Unit = field(init=False)
    temperature       : Unit = field(init=False)
    pressure          : Unit = field(init=False)
    dynamic_viscosity : Unit = field(init=False)
    charge            : Unit = field(init=False)
    dipole_moment     : Unit = field(init=False)
    electric_field    : Unit = field(init=False)
    density           : Unit = field(init=False)

    dt            : Quantity = field(init=False)
    skin          : Quantity = field(init=False)

# Concrete unit style reference implementations
class LAMMPSUnitStyleReal(LAMMPSUnitStyle):
    STYLE_NAME = 'real'
    __doc__ = f'''Unit preferences for the "{STYLE_NAME}" unit style''' # TODO : find way for this to be populated during inheritance

    mass              = gram / mole
    distance          = angstrom
    time              = femtosecond
    energy            = kilocalorie_per_mole
    velocity          = angstrom / femtosecond
    force             = kilocalorie_per_mole / angstrom
    torque            = kilocalorie_per_mole
    temperature       = kelvin
    pressure          = atmosphere
    dynamic_viscosity = poise
    charge            = elementary_charge # same as electron charge
    dipole_moment     = elementary_charge * angstrom
    electric_field    = volt / angstrom
    density           = gram / centimeter**3 # the LAMMPS docs have the length exponent as "dim" (how to make this dimension-dependent?)

    dt   = 1.0 * femtosecond
    skin = 2.0 * angstrom

class LAMMPSUnitStyleMetal(LAMMPSUnitStyle):
    STYLE_NAME = 'metal'
    __doc__ = f'''Unit preferences for the "{STYLE_NAME}" unit style'''

    mass              = gram / mole
    distance          = angstrom
    time              = picosecond
    energy            = electronvolt
    velocity          = angstrom / picosecond
    force             = electronvolt / angstrom
    torque            = electronvolt
    temperature       = kelvin
    pressure          = bar
    dynamic_viscosity = poise
    charge            = elementary_charge
    dipole_moment     = elementary_charge * angstrom
    electric_field    = volt / angstrom
    density           = gram / centimeter**3 

    dt   = 1E-3 * picosecond
    skin = 2.0 * angstrom

class LAMMPSUnitStyleSI(LAMMPSUnitStyle):
    STYLE_NAME = 'si'
    __doc__ = f'''Unit preferences for the "{STYLE_NAME}" unit style'''

    mass              = kilogram
    distance          = meter
    time              = second
    energy            = joule
    velocity          = meter / second
    force             = newton
    torque            = newton * meter
    temperature       = kelvin
    pressure          = pascal
    dynamic_viscosity = pascal * second
    charge            = coulomb
    dipole_moment     = coulomb * meter
    electric_field    = volt / meter
    density           = kilogram * meter**3

    dt   = 1E-8 * second
    skin = 1E-3 * meter

class LAMMPSUnitStyleCGS(LAMMPSUnitStyle):
    STYLE_NAME = 'cgs'
    __doc__ = f'''Unit preferences for the "{STYLE_NAME}" unit style'''

    mass              = gram
    distance          = centimeter
    time              = second
    energy            = erg
    velocity          = centimeter / second
    force             = dyne
    torque            = dyne * centimeter
    temperature       = kelvin
    pressure          = dyne / centimeter**2
    dynamic_viscosity = poise
    charge            = statcoulomb
    dipole_moment     = statcoulomb * centimeter
    electric_field    = statvolt / centimeter
    density           = gram / centimeter**3 

    dt   = 1E-8 * second
    skin = 0.1 * centimeter

class LAMMPSUnitStyleElectron(LAMMPSUnitStyle):
    STYLE_NAME = 'electron'
    __doc__ = f'''Unit preferences for the "{STYLE_NAME}" unit style'''
    
    mass = dalton
    distance = bohr
    time = femtosecond
    energy = hartree
    velocity = bohr / atomic_time # LAMMPS docs definition of atomic time seems to differ from convention (given as 1.03275e-15 seconds by LAMMPS)
    force = hartree / bohr
    temperature = kelvin
    pressure = pascal
    charge = elementary_charge
    dipole_moment = debye # = statcoulomb * centimeter # equivalent to one Debye
    electric_field = volt / centimeter
    density = None # this is not defined in the LAMMPS docs (why??)

    dt   = 1E-3 * femtosecond
    skin = 2.0 * bohr

class LAMMPSUnitStyleMicro(LAMMPSUnitStyle):
    STYLE_NAME = 'micro'
    __doc__ = f'''Unit preferences for the "{STYLE_NAME}" unit style'''

    mass = picogram
    distance = micrometer
    time = microsecond
    energy = picogram * micrometer**2 / microsecond**2
    velocity = micrometer / microsecond
    force = picogram * micrometer / microsecond**2
    torque = picogram * micrometer**2 / microsecond**2
    temperature = kelvin
    pressure = picogram / (micrometer * microsecond**2)
    dynamic_viscosity = picogram / (micrometer * microsecond)
    charge = picocoulomb
    dipole = picocoulomb * micrometer
    electric_field = volt / micrometer
    density = picogram / micrometer**3

    dt   = 2.0 * microsecond
    skin = 0.1 * micrometer

class LAMMPSUnitStyleNano(LAMMPSUnitStyle):
    STYLE_NAME = 'nano'
    __doc__ = f'''Unit preferences for the "{STYLE_NAME}" unit style'''

    mass              = attogram
    distance          = nanometer
    time              = nanosecond
    energy            = (attogram * nanometer**2) / nanosecond**2
    velocity          = nanometer / nanosecond
    force             = (attogram * nanometer) / nanosecond**2
    torque            = (attogram * nanometer**2) / nanosecond**2
    temperature       = kelvin
    pressure          = attogram / (nanometer * nanosecond**2)
    dynamic_viscosity = attogram / (nanometer * nanosecond)
    charge            = elementary_charge
    dipole            = charge * nanometer
    electric_field    = volt / nanometer
    density           = attogram / nanometer**3

    dt   = 4.5E-4 * nanosecond
    skin = 0.1 * nanometer

# TODO : implement LJ-style? (is mostly dimensionless, with dt = 0.005 tau and skin = 0.3 sigma)