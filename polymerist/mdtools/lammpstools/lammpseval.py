'''For extracting information from LAMMPS input and data files'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Optional

import lammps
from openmm.unit import Unit, Quantity
from openmm.unit import degree, angstrom, kilojoule_per_mole

from .unitstyles import LAMMPSUnitStyle
from ...genutils.decorators.functional import allow_pathlib_paths


# CONSTANTS
SUPPRESS_LAMMPS_STDOUT : list[str] = ["-screen", "none", "-log", "none"]
LAMMPS_ENERGY_KW_ALIASES : dict[str, str] = {
    # bonded
    'E_bond'    : 'Bond',
    'E_angle'   : 'Angle',
    'E_dihed'   : 'Proper Torsion',
    'E_impro'   : 'Improper Torsion',
    'E_mol'     : 'Bonded', # = bond + angle + torsion + impropers
    # nonbonded
    'E_coul'    : 'Coulomb Short',
    'E_long'    : 'Coulomb Long',
    'E_vdwl'    : 'vdW',
    'E_tail'    : 'Dispersion',
    'E_pair'    : 'Nonbonded',
    # overall
    'PotEng'    : 'Potential',
    'KinEng'    : 'Kinetic',
    'TotEng'    : 'Total',  # = KE + PE
    'Enthalpy'  : 'Enthalpy', # = total + PV work
    'Ecouple'   : 'Thermostat', # losses from thermostat/barostat action
    'Econserve' : 'Unconserved',  # deviation from ideality, = total + thermostat
}

# READER FUNCTIONS
@allow_pathlib_paths
def get_lammps_unit_style(
        lmp_input_path : str,
        cmdargs : Optional[list[str]]=None,
    ) -> LAMMPSUnitStyle:
    '''Fetch Pythonic representation of unit styes specified for a LAMMPS input file'''
    cmdargs = [] if cmdargs is None else [arg for arg in cmdargs] # make copy to shield input list from being mangled by lammps.lammps() (adds extra stuff to this list)
    with lammps.lammps(cmdargs=cmdargs) as lmp:
        lmp.file(lmp_input_path)
        unit_style_name : str = lmp.extract_global('units')
        
    try:
        return LAMMPSUnitStyle.subclass_registry[unit_style_name]
    except KeyError:
        raise ValueError(f'Invalid unit style "{unit_style_name}" specified in LAMMPS input file')
    
@allow_pathlib_paths
def get_lammps_lattice_parameters(
        lmp_input_path : str,
        cmdargs : Optional[list]=None,
        preferred_length_unit : Unit=None,
        preferred_angle_unit : Unit=None,
    ) -> dict[str, Quantity]: # TODO: reimplement this with maths.lattices.bravais.LatticeParameters
    '''Extract the 6 lattice parameters (i.e. 3 lattice vector lengths and 3 inter-vector angles)
    specified for the simulation box defined by a LAMMPS input file'''
    # sanitize unit preferences
    ## lengths
    calculated_length_unit : Unit = get_lammps_unit_style(lmp_input_path, cmdargs=cmdargs).distance
    if preferred_length_unit is None:
        preferred_length_unit = calculated_length_unit
    assert(preferred_length_unit.is_compatible(angstrom))
    
    ## angles
    calculated_angle_unit : Unit = degree # DEVNOTE: AFAIK, LAMMPS outputs lattice vector angles in degrees (thought couldn't find documentation explicitly stating this)
    if preferred_angle_unit is None:
        preferred_angle_unit = calculated_angle_unit
    assert(preferred_angle_unit.is_compatible(degree))
    
    lattice_param_units : dict[str, tuple[Unit, Unit]] = { # indicates unit LAMMPS will return, and the unit preferred
        'cella' : (calculated_length_unit, preferred_length_unit),
        'cellb' : (calculated_length_unit, preferred_length_unit),
        'cellc' : (calculated_length_unit, preferred_length_unit),
        'cellalpha' : (calculated_angle_unit, preferred_angle_unit),
        'cellbeta'  : (calculated_angle_unit, preferred_angle_unit),
        'cellgamma' : (calculated_angle_unit, preferred_angle_unit),
    }
    
    # read lattice parameters from input
    cmdargs = [] if cmdargs is None else [arg for arg in cmdargs] # make copy to shield input list from being mangled by lammps.lammps() (adds extra stuff to this list)
    with lammps.lammps(cmdargs=cmdargs) as lmp:
        lmp.file(lmp_input_path)
        return {
            thermo_kw : (lmp.get_thermo(thermo_kw)*calculated_unit).in_units_of(preferred_unit)
                for thermo_kw, (calculated_unit, preferred_unit) in lattice_param_units.items()
        }
        
@allow_pathlib_paths
def get_lammps_energies(
        lmp_input_path : str,
        cmdargs : Optional[list]=None,
        preferred_energy_unit : Optional[Unit]=None,
        energy_kw_remap : Optional[dict[str, str]]=None,
    ) -> dict[str, Quantity]:
    '''Perform single-point energy evaluation from a LAMMPS input file, respecting the LAMMPS unit style specified in the input'''
    if energy_kw_remap is None:
        energy_kw_remap = LAMMPS_ENERGY_KW_ALIASES
    
    # sanitize unit preferences
    calculated_energy_unit : Unit = get_lammps_unit_style(lmp_input_path, cmdargs=cmdargs).energy
    if preferred_energy_unit is None:
        preferred_energy_unit = calculated_energy_unit
    assert(preferred_energy_unit.is_compatible(kilojoule_per_mole))
    
    # read lattice parameters from input
    cmdargs = [] if cmdargs is None else [arg for arg in cmdargs] # make copy to shield input list from being mangled by lammps.lammps() (adds extra stuff to this list)
    with lammps.lammps(cmdargs=cmdargs) as lmp:
        lmp.file(lmp_input_path)
        # lmp.command('run 0') # ensure single-step evaluation is executed
        
        return {
            energy_kw_remap.get(energy_kw, energy_kw) : (energy_magnitude*calculated_energy_unit).in_units_of(preferred_energy_unit)
                for energy_kw, energy_magnitude in lmp.last_thermo().items()
        }