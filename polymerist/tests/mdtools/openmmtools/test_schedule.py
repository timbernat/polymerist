'''Behavioral tests for OpenMM serial simulation schedule runner'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import pytest

from numpy.typing import NDArray

from openmm.unit import (
    Quantity as OpenMMQuantity,
    femtosecond, picosecond,
    nanometer, centimeter,
    gram, kelvin, atmosphere,
)
from openmm import System, State
from openmm.app import Topology as OpenMMTopology

from openff.toolkit import Molecule, ForceField
from openff.interchange import Interchange

from polymerist.mdtools.openfftools.boxvectors import get_topology_bbox, pad_box_vectors_uniform
from polymerist.mdtools.openfftools.solvation.packing import pack_topology_with_solvent
from polymerist.mdtools.openfftools.solvation.solvents import water_TIP3P

from polymerist.mdtools.openmmtools.parameters import (
    ThermoParameters,
    ThermostatParameters,
    BarostatParameters,
    IntegratorParameters,
    ReporterParameters,
    SimulationParameters,
)


# FIXTURES
@pytest.fixture(scope='module')
def interchange() -> Interchange:
    '''OpenFF Interchange object from which to derive OpenMM simulation inputs'''
    offmol = Molecule.from_smiles('c1ccccc1C(=O)O')
    offmol.generate_conformers(n_conformers=1)
    offmol.assign_partial_charges(partial_charge_method='gasteiger')
    offtop = offmol.to_topology()

    # NOTE: need to pack box to ensure box does not collapse below nonbonded cutoff during NPT simulation tests 
    bbox = pad_box_vectors_uniform(
        get_topology_bbox(offtop),
        pad_amount=0.9*nanometer,
    ) # guarantee box is bigger than 1 nm cutoff in all directions
    packed_top = pack_topology_with_solvent(
        offtop,
        solvent=water_TIP3P,
        box_vecs=bbox,
        density=1.0*(gram*centimeter**-3),
    )
    
    ff = ForceField('openff-2.2.0.offxml')
    inc = ff.create_interchange(packed_top, charge_from_molecules=[offmol, water_TIP3P])
    inc.box = bbox
    
@pytest.fixture(scope='function') # DEVNOTE: should be re-made per call to avoid coss-contamination
def openmm_topology(interchange : Interchange) -> OpenMMTopology:
    '''OpenMM Topology object for the simulation'''
    return interchange.to_openmm_topology(collate=False)

@pytest.fixture(scope='function') # DEVNOTE: should be re-made per call to avoid coss-contamination
def openmm_system(interchange : Interchange) -> System:
    '''OpenMM System object for the simulation'''
    return interchange.to_openmm_system(combine_nonbonded_forces=False)

@pytest.fixture(scope='function') # DEVNOTE: should be re-made per call to avoid coss-contamination
def openmm_positions(interchange : Interchange) -> OpenMMQuantity:
    '''OpenMM positions for the simulation'''
    return interchange.get_positions(include_virtual_sites=True).to_openmm()

## TODO: load optional compatible initial state to test initial State injection

@pytest.fixture(scope='module')
def sim_params_NPT() -> SimulationParameters:
    '''Simulation parameters for thermostatted Simulation'''
    return SimulationParameters(
        thermo_params=ThermoParameters(
            thermostat_params=ThermostatParameters(
                temperature=500*kelvin,
                timescale=1.5/picosecond,
                thermostat='LangevinMiddle',
            ),
            barostat_params=BarostatParameters(
                pressure=1.0*atmosphere,
                temperature=500,
                update_frequency=25,
                barostat='MonteCarlo',
            ),
        ),
        integ_params=IntegratorParameters(
            time_step=1.0*femtosecond,
            total_time=10*picosecond,
            num_samples=10,
        ),
        reporter_params=ReporterParameters(
            report_checkpoint=True,
            report_state=True,
            report_trajectory=True,
            report_state_data=True,
        ),
    )
    
@pytest.fixture(scope='module')
def sim_params_NVT() -> SimulationParameters:
    '''Simulation parameters for thermostatted AND barostatted Simulation'''
    return SimulationParameters(
        thermo_params=ThermoParameters(
            thermostat_params=ThermostatParameters(
                temperature=300*kelvin,
                timescale=1.0/picosecond,
                thermostat='Andersen',
            ),
        ),
        integ_params=IntegratorParameters(
            time_step=1.0*femtosecond,
            total_time=5*picosecond,
            num_samples=50,
        ),
        reporter_params=ReporterParameters(
            report_checkpoint=True,
            report_state=True,
            report_trajectory=True,
            report_state_data=True,
        ),
    )
# TODO: add premade State compatible with the System defined here to test state injection at start of schedule
    
# TESTS PROPER