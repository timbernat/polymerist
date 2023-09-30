'''Boilerplate for setting up OpenMM Simulations and related files'''

import logging
LOGGER = logging.getLogger(__name__)

from typing import Iterable, Optional, Union
from pathlib import Path

from openmm import Force, System, State, XmlSerializer
from openmm.app import Simulation, Topology
from openmm.unit import Quantity

# from .records import SimulationParameters
from .serialization import serialize_topology_from_simulation, serialize_system, SimulationPaths
from .thermo import ThermoParameters, EnsembleFactory
from .reporters import ReporterParameters
from .parameters import IntegratorParameters


# TODO : add optional_in_place??
def label_forces(system : System) -> None:
    '''Designates each Force in a System with a unique force group and assigns helpful names by Force type'''
    for i, force in enumerate(system.getForces()):
        force.setForceGroup(i)

    # TODO : add labelling (depends partially on Interchange's NonbondedForce separation)

def simulation_from_thermo(topology : Topology, system : System, thermo_params : ThermoParameters, time_step : Quantity, state : Optional[State]=None) -> Simulation:
    '''Prepare an OpenMM simulation from a serialized thermodynamics parameter set'''
    ens_fac = EnsembleFactory.from_thermo_params(thermo_params)
    if (forces := ens_fac.forces()):
        for force in forces:
            system.addForce(force) # add forces to System BEFORE creating Simulation to avoid having to reinitialze the Conext to preserve changes 
    label_forces(system) # ensure all system forces (including any ensemble-specific ones) are labelled

    simulation = Simulation(
        topology=topology,
        system=system,
        integrator=ens_fac.integrator(time_step),
        state=state
    )

    return simulation

def initialize_simulation_and_files(out_dir : Path, out_name : str, sim_paths : SimulationPaths, topology : Topology, system : System) -> Simulation:
    '''Create simulation, bind Reporters, and update simulation Paths with newly-generated files'''
    # extract info from 
    integ_params = IntegratorParameters.from_file( sim_paths.integ_params)
    thermo_params = ThermoParameters.from_file(    sim_paths.thermo_params)
    reporter_params = ReporterParameters.from_file(sim_paths.reporter_params)

    # if sim_paths.state is None:
    #     state = None
    # else:
    try:
        with sim_paths.state.open('r') as state_file:
            state = XmlSerializer.deserialize(state_file.read())
    except (AttributeError, ValueError): # handle cases where state path is None or state file is empty, respectively
        state = None

    # create simulation and add reporters
    simulation = simulation_from_thermo(topology, system, thermo_params, time_step=integ_params.time_step, state=state)
    rep_paths, reps = reporter_params.prepare_reporters(out_dir, out_name, integ_params.report_interval)
    
    for reporter in reps:
        simulation.reporters.append(reporter) # add reporters to simulation instance

    for attr, path in rep_paths.items():
        setattr(sim_paths, attr, path) # register reporter output paths

    return simulation


def record_simulation_top_and_sys(out_dir : Path, out_name : str, simulation : Simulation, sim_paths : SimulationPaths) -> None:
    '''Serilaize current topology and initial state, system, and checkpoint'''
    sim_paths.topology = serialize_topology_from_simulation(simulation, out_dir, out_name) # TODO : make this consistent with rest of Path output
    LOGGER.info(f'Saved energy-minimized Simulation Topology at {sim_paths.topology}')

    sim_paths.system = serialize_system(simulation.system, out_dir, out_name)
    LOGGER.info(f'Saved serialized Simulation System at {sim_paths.system}')
