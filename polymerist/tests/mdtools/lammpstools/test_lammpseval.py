'''Unit tests for `lammpseval` package'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import pytest
from pathlib import Path

import polymerist.tests.data as testdata
from polymerist.genutils.fileutils.filetree import temporary_cd
from polymerist.genutils.importutils.pkginspect import get_dir_path_within_package

from openmm.unit import kilojoule_per_mole, kilocalorie_per_mole, Quantity, angstrom, degree
from polymerist.mdtools.lammpstools.unitstyles import LAMMPSUnitStyleReal
from polymerist.mdtools.lammpstools.lammpseval import (
    get_lammps_unit_style,
    get_lammps_lattice_parameters,
    get_lammps_energies,
    SUPPRESS_LAMMPS_STDOUT,
)


@pytest.fixture
def lammps_dir() -> Path:
    return get_dir_path_within_package('LAMMPS_files', testdata)

@pytest.fixture
def lammps_input_path() -> Path:
    return Path('inputs.in') # important to note that this Path is RELATIVE to the LAMMPS file directory

@pytest.fixture
def lammps_data_path() -> Path:
    return Path('data.in') # important to note that this Path is RELATIVE to the LAMMPS file directory


def test_get_unit_style(lammps_dir : Path, lammps_input_path : Path) -> None:
    '''Test that correct unit style object is read from LAMMPS files'''
    with temporary_cd(lammps_dir): # necessary to get LAMMPS to relate the input and data correctly
        assert get_lammps_unit_style(lammps_input_path, cmdargs=SUPPRESS_LAMMPS_STDOUT) == LAMMPSUnitStyleReal
        
def test_get_lattice_parameters(lammps_dir : Path, lammps_input_path : Path) -> None:
    '''Test that correct lattice parameters are read from LAMMPS files'''
    with temporary_cd(lammps_dir): # necessary to get LAMMPS to relate the input and data correctly
        assert get_lammps_lattice_parameters(
            lammps_input_path,
            cmdargs=SUPPRESS_LAMMPS_STDOUT,
            preferred_length_unit=angstrom,
            preferred_angle_unit=degree,
        ) == {
            'cella' : 12.357899665*angstrom,
            'cellb' : 15.890199661*angstrom,
            'cellc' : 13.925300598*angstrom,
            'cellalpha' : 90.0*degree,
            'cellbeta'  : 90.0*degree,
            'cellgamma' : 90.0*degree,
        }
        
def test_get_energies(lammps_dir : Path, lammps_input_path : Path) -> None:
    '''Test that single-point energy evaluation is calculated and labelled correctly'''
    with temporary_cd(lammps_dir): # necessary to get LAMMPS to relate the input and data correctly
        assert get_lammps_energies(
            lammps_input_path,
            cmdargs=SUPPRESS_LAMMPS_STDOUT,
            preferred_energy_unit=kilocalorie_per_mole,
        ) == {
            'Bond': 29.170789900387238*kilocalorie_per_mole,
            'Angle': 330.9716442076586*kilocalorie_per_mole,
            'Proper Torsion': 79.17111534664818*kilocalorie_per_mole,
            'Improper Torsion': 0.18343504929981924*kilocalorie_per_mole,
            'Nonbonded': 433093.4641287483*kilocalorie_per_mole,
            'vdW': 433528.89352945634*kilocalorie_per_mole,
            'Coulomb Short': 340.6987993890613*kilocalorie_per_mole,
            'Coulomb Long': -776.128200097151*kilocalorie_per_mole,
            'Dispersion': -4.779493985001752*kilocalorie_per_mole,
            'Potential': 433532.9611132523*kilocalorie_per_mole,
        }
        