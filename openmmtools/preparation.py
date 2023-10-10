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
from .thermo import EnsembleFactory
from .parameters import ThermoParameters, ReporterParameters, IntegratorParameters, SimulationParameters


def label_forces(system : System) -> None:
    '''Designates each Force in a System with a unique force group and assigns helpful names by Force type'''
    for i, force in enumerate(system.getForces()):
        force.setForceGroup(i)

    # TODO : add labelling (depends partially on Interchange's NonbondedForce separation)

def simulation_from_thermo(topology : Topology, system : System, thermo_params : ThermoParameters, time_step : Quantity, state : Optional[Path]=None) -> Simulation:
    '''Prepare an OpenMM simulation from a serialized thermodynamics parameter set'''
    ens_fac = EnsembleFactory.from_thermo_params(thermo_params)
    if (forces := ens_fac.forces()): # check if any extra forces are present
        for force in forces:
            system.addForce(force) # add forces to System BEFORE creating Simulation to avoid having to reinitialize the Conext to preserve changes 
            LOGGER.info(f'Added {force.getName()} Force to System')
    label_forces(system) # ensure all system forces (including any ensemble-specific ones) are labelled

    if state is not None:
        try:
            with state.open('r') as statefile:
                saved_state = XmlSerializer.deserialize(statefile.read())
        except ValueError: # catch when a state file exists but is invalid (or indeed empty)
            state = None

    simulation = Simulation(
        topology=topology,
        system=system,
        integrator=ens_fac.integrator(time_step),
        state=state
    )

    return simulation

def initialize_simulation_and_files(out_dir : Path, prefix : str, sim_params : SimulationParameters, topology : Topology, system : System, positions : Optional[Quantity]=None) -> tuple[Simulation, SimulationPaths]:
    '''Create simulation, bind Reporters, and update simulation Paths with newly-generated files'''
    sim_paths = SimulationPaths.from_dir_and_parameters(out_dir, prefix, sim_params, touch=True)

    # create simulation and add reporters
    simulation = simulation_from_thermo(topology, system, sim_params.thermo_params, time_step=sim_params.integ_params.time_step, state=sim_paths.state_path)
    for reporter in sim_params.reporter_params.prepare_reporters(report_interval=sim_params.integ_params.report_interval):
        simulation.reporters.append(reporter) # add reporters to simulation instance

    if positions is not None:
        # TODO : add position type/shape checking
        LOGGER.info('Setting positions in context')
        simulation.context.setPositions(positions) # set positions if provided

    return simulation, sim_paths

def run_simulation_schedule(working_dir : Path, schedule : dict[str, SimulationParameters], init_top : Topology, init_sys : System, init_pos : Quantity) -> None:
    '''Run several OpenMM simulations in series, based on an initial set of OpenMM objects and a "schedule" consisting of a sequence of named parameter sets'''
    working_dir.mkdir(exist_ok=True)

    num_steps = len(schedule)
    for i, (step_name, sim_params) in enumerate(schedule.items()): # TOSELF : may want to shift to an explicitly-ordered map (e.g. collections OrderedDict); Python 3.6+ preserves order, but wouldn't hurt to be extra safe
        if i == 0:
            ommtop, ommsys, ommpos = init_top, init_sys, init_pos # use initial Topology and System for first sim
        else:
            ommtop, ommsys, ommpos = ommsim.topology, ommsim.system, ommsim.context.getState(getPositions=True).getPositions(asNumpy=True) # use Topology and System from previous sim for next sim

        LOGGER.info(f'Initializing simulation {i + 1}/{num_steps} ({step_name})')
        ommsim, sim_paths = initialize_simulation_and_files(
            out_dir=working_dir / step_name,
            prefix=step_name,
            sim_params=sim_params,
            topology=ommtop,
            system=ommsys,
            positions=ommpos
        )
        
        LOGGER.info('Performing energy minimization')
        ommsim.minimizeEnergy()
        LOGGER.info('Energy successfully minimized')

        serialize_topology_from_simulation(sim_paths.topology_path, ommsim) # TODO : make this consistent with rest of Path output
        LOGGER.info(f'Saved energy-minimized Simulation Topology at {sim_paths.topology_path}')

        serialize_system(sim_paths.system_path, ommsim.system)
        LOGGER.info(f'Saved serialized Simulation System at {sim_paths.system_path}')

        LOGGER.info(f'Integrating {sim_params.integ_params.total_time} OpenMM Simulation for {sim_params.integ_params.num_steps} steps')
        ommsim.step(sim_params.integ_params.num_steps)
        LOGGER.info('Simulation integration completed successfully')