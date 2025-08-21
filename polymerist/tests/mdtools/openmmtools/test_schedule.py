'''Behavioral tests for OpenMM serial simulation schedule runner'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import pytest
from typing import Union

import numpy as np
import pandas as pd

from openmm.unit import (
    Quantity as OpenMMQuantity,
    femtosecond, picosecond,
    nanometer, centimeter,
    gram, kelvin, atmosphere,
)
from openmm import System, State, Context
from openmm.app import Simulation, Topology as OpenMMTopology

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
from polymerist.mdtools.openmmtools.serialization.paths import SimulationPaths
from polymerist.mdtools.openmmtools.serialization.state import DEFAULT_STATE_PROPS
from polymerist.mdtools.openmmtools.execution import run_simulation_schedule


# FIXTURES
@pytest.fixture(scope='session')
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
    
    return inc
    
 # DEV: would have love to make this per-function to avoid cross-contamination,
 # but the session-scoped simulations below force me to make these session wide as well
 # (will raise ScopeError and not run test if I don't do this)
@pytest.fixture(scope='session')
def openmm_topology(interchange : Interchange) -> OpenMMTopology:
    '''OpenMM Topology object for the simulation'''
    return interchange.to_openmm_topology(collate=False)

@pytest.fixture(scope='session')
def openmm_system(interchange : Interchange) -> System:
    '''OpenMM System object for the simulation'''
    return interchange.to_openmm_system(combine_nonbonded_forces=False)

@pytest.fixture(scope='session')
def openmm_positions(interchange : Interchange) -> OpenMMQuantity:
    '''OpenMM positions for the simulation'''
    return interchange.get_positions(include_virtual_sites=True).to_openmm()

@pytest.fixture(scope='session')
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
            total_time=0.5*picosecond, # NOTE: deliberately running super short - not testing quality of sims, just how they run
            num_samples=10,
        ),
        reporter_params=ReporterParameters(
            report_checkpoint=True,
            report_state=True,
            report_trajectory=True,
            report_state_data=True,
        ),
    )
    
@pytest.fixture(scope='session')
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
            total_time=0.4*picosecond, # NOTE: deliberately running super short - not testing quality of sims, just how they run
            num_samples=20,
        ),
        reporter_params=ReporterParameters(
            report_checkpoint=True,
            report_state=True,
            report_trajectory=True,
            report_state_data=True,
        ),
    )
    
@pytest.fixture(scope='session') # NOTE: running simulations is pretty expensive and should only be done once
def simulation_history(
        tmp_path_factory, # this is a "magic" keyword from PyTest to generate a session-long temporary Path
        sim_params_NPT : SimulationParameters,
        sim_params_NVT : SimulationParameters,
        openmm_topology : OpenMMTopology,
        openmm_system : System,
        openmm_positions : np.ndarray[float],
    ) -> dict[str, dict[str, Union[Simulation, SimulationPaths]]]: 
    '''The outputs of short simulations specified by the parameters above run in order'''
    return run_simulation_schedule(
        working_dir=tmp_path_factory.mktemp('openmm'),
        schedule={
            'NPT' : sim_params_NPT,
            'NVT' : sim_params_NVT,
        },
        init_top=openmm_topology,
        init_sys=openmm_system,
        init_pos=openmm_positions,
        init_state=None, # DEV: for now
        return_history=True
    )


# TESTS
def test_barostat_NPT(simulation_history) -> None:
    '''Test that barostat action is being applied to constant-pressure simulations'''
    sim_paths = simulation_history['NPT']['paths']
    state_data : pd.DataFrame = pd.read_csv(sim_paths.state_data_path)
    box_volumes : np.ndarray[float] = state_data['Box Volume (nm^3)'].values

    assert not np.allclose(box_volumes[0], box_volumes) # simulation volume SHOULD be changing (i.e. non-constant) to maintain pressure

# DEV: the way the simulation schedule is set up here ALSO tests that no barostat bleedover occurs (NPT is followed by NVT)
def test_barostat_NVT_no_bleedover(simulation_history) -> None:
    '''Test that barostat action is being applied to constant-pressure simulations'''
    sim_paths = simulation_history['NVT']['paths']
    state_data : pd.DataFrame = pd.read_csv(sim_paths.state_data_path)
    box_volumes : np.ndarray[float] = state_data['Box Volume (nm^3)'].values

    assert np.allclose(box_volumes[0], box_volumes) # simulation volume should be constant

def test_inject_state(
    tmp_path_factory, # this is a "magic" keyword from PyTest to generate a session-long temporary Path
    sim_params_NVT : SimulationParameters,
    openmm_topology : OpenMMTopology,
    openmm_system : System,
    openmm_positions : np.ndarray[float],
) -> None:
    '''Test that injection of custom state into first Simulation in schedule works'''
    # pre-load initial State to use for Simulation - DEV: loading state from file proved too unreliable for CI tests
    context = Context(
        openmm_system,
        sim_params_NVT.thermo_params.integrator(time_step=sim_params_NVT.integ_params.time_step),
    )
    context.setPositions(openmm_positions)
    state = context.getState(**DEFAULT_STATE_PROPS) 
    
    run_simulation_schedule(
        working_dir=tmp_path_factory.mktemp('openmm'),
        schedule={'sim' : sim_params_NVT},
        init_top=openmm_topology,
        init_sys=openmm_system,
        init_pos=openmm_positions,
        init_state=state,
        return_history=True
    )
    # no explicit assert - just making sure no Exceptions are raised
