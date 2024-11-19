'''For gathering information and running calculations from LAMMPS input files'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Optional, Union

import re
from pathlib import Path

import lammps
from openmm.unit import Unit, Quantity
from openmm.unit import degree, kilocalorie_per_mole

from .unitstyles import LAMMPSUnitStyle
from ...genutils.decorators.functional import allow_string_paths, allow_pathlib_paths


_DISTANCE_WILDCARD : str = '$DISTANCE'
LMP_UNIT_REGEX = re.compile(r'^units\s(?P<unit_style>\w*)$')
LMP_THERMO_STYLE_REGEX = re.compile(r'^thermo_style\s(?P<thermo_style>\b\w*?\b)\s(?P<energy_evals>.*)$')

LAMMPS_ENERGY_KW : dict[str, str] = { # keywordS for evaluating energies and their corresponding meanings
    'ebond'  : 'Bond',
    'eangle' : 'Angle',
    'edihed' : 'Proper Torsion',
    'eimp'   : 'Improper Torsion',
    'ecoul'  : 'Coulomb Short',
    'elong'  : 'Coulomb Long',
    'evdwl'  : 'vdW',
    'etail'  : 'Dispersion',
    'epair'  : 'Nonbonded',
    'pe'     : 'Potential',
    'ke'     : 'Kinetic',
    'etotal' : 'Total'
}

LAMMPS_CELL_KW = ( # keywords for probing unit cell sizes and angles
    'cella', # TODO : add aliases
    'cellb',
    'cellc',
    'cellalpha',
    'cellbeta',
    'cellgamma',
)

LAMMPS_CELL_UNITS : dict[str, Union[str, Unit]] = { # units associated with unit cell keywords
    'cella' : _DISTANCE_WILDCARD, 
    'cellb' : _DISTANCE_WILDCARD,
    'cellc' : _DISTANCE_WILDCARD,
    'cellalpha' : degree,
    'cellbeta'  : degree,
    'cellgamma' : degree,
}

@allow_string_paths
def parse_lammps_input(lmp_in_path : Path) -> dict[str, Union[str, LAMMPSUnitStyle]]: # NOTE : this can and will be expanded in the future
    '''Read unit style, thermo style, and energies to evaluate from a LAMMPS input file'''
    info_dict = {}
    with lmp_in_path.open('r') as lmp_in_file:
        for line in lmp_in_file.read().split('\n'):
            if (thermo_match := re.search(LMP_THERMO_STYLE_REGEX, line)):
                info_dict.update(thermo_match.groupdict())
                info_dict['energy_evals'] = info_dict['energy_evals'].split(' ') # separate on spaces (TODO : maybe find more elegant way to do this in the future?)
            
            if (units_match := re.search(LMP_UNIT_REGEX, line)):
                unit_style_name = units_match.groupdict()['unit_style']
                try:
                    info_dict['unit_style'] = LAMMPSUnitStyle.subclass_registry[unit_style_name]
                except KeyError:
                    raise ValueError(f'Invalid unit style "{unit_style_name}" specified in LAMMPS input file')
                
    return info_dict

@allow_pathlib_paths
def get_lammps_energies(lmp_in_path : str, preferred_unit : Unit=kilocalorie_per_mole, cmdargs : Optional[list]=None) -> dict[str, Quantity]:
    '''Perform an energy evaluation using a LAMMPS input file
    Alternative to interchange.drivers.get_lammps_energies which is dynamically aware of energy units and assumes nothing about which energies are specified by thermo_style'''
    cmdargs = [] if cmdargs is None else [arg for arg in cmdargs] # need to copy list to keep read-only (lammps.lammps() modifies the argument list passed in-place)
    
    assert(preferred_unit.is_compatible(kilocalorie_per_mole)) # whatever unit is desired, it must be one of energy
    lammps_info = parse_lammps_input(lmp_in_path)
    energy_unit     = lammps_info['unit_style'].energy
    energy_contribs = lammps_info['energy_evals']

    with lammps.lammps(cmdargs=cmdargs) as lmp: # need to create new lammps() object instance for each run
        # lmp.commands_string( ENERGY_EVAL_STR.replace('$INP_FILE', str(lammps_file)) )
        lmp.file(lmp_in_path) # read input file and calculate energies; NOTE that this NEEDS to be a string (not Path!)
        return {
            f'{LAMMPS_ENERGY_KW[contrib]}' : (lmp.get_thermo(contrib) * energy_unit).in_units_of(preferred_unit)
                for contrib in energy_contribs
        }

@allow_pathlib_paths
def get_lammps_unit_cell(lmp_in_path : str, cmdargs : Optional[list]=None) -> dict[str, Quantity]: # TODO: reimplement this with maths.lattices.bravais.LatticeParameters
    '''Extract the 6 unit cell parameters specified by a LAMMPS input file'''
    cmdargs = [] if cmdargs is None else [arg for arg in cmdargs] # need to copy list to keep read-only (lammps.lammps() modifies the argument list passed in-place)
    
    lammps_info = parse_lammps_input(lmp_in_path)
    length_unit = lammps_info['unit_style'].distance

    with lammps.lammps(cmdargs=cmdargs) as lmp: # need to create new lammps() object instance for each run
        lmp.file(lmp_in_path) # read input file and calculate energies; NOTE that this NEEDS to be a string (not Path!)
        cell_params = {}
        for cell_param_kw in LAMMPS_CELL_KW:
            cell_param_unit = LAMMPS_CELL_UNITS[cell_param_kw]
            if cell_param_unit == _DISTANCE_WILDCARD:
                cell_param_unit = length_unit

            cell_params[cell_param_kw] = lmp.get_thermo(cell_param_kw) * cell_param_unit
    
    return cell_params