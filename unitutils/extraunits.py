'''Defining units which, for one reason or another, are not defined in Pint or OpenMM units'''

from openmm.unit import Unit, ScaledUnit 

# DEFINING ELECTRON VOLTS
from openmm.unit import joule
from scipy.constants import electron_volt as electron_volt_joules

electronvolt_base = ScaledUnit(electron_volt_joules, joule, 'electronvolt', 'eV')
electronvolt = eV = Unit({electronvolt_base : 1.0})
