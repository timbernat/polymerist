'''Boilerplate for setting up OpenMM Simulations and related files'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import logging
LOGGER = logging.getLogger(__name__)

from typing import Optional
from pathlib import Path

from openmm import System
from openmm.app import Simulation, Topology
from openmm.unit import Quantity

from .thermo import ThermoParameters
from .parameters import SimulationParameters
from .serialization.paths import SimulationPaths
from .serialization.state import StateLike, load_state_flexible
from .forces import _POLYMERIST_FORCE_GROUP, force_added_by_polymerist, impose_unique_force_groups


def simulation_from_thermo(
        topology : Topology,
        system : System,
        thermo_params : ThermoParameters,
        time_step : Quantity,
        positions : Optional[Quantity]=None,
        state : Optional[StateLike]=None,
    ) -> Simulation:
    '''Prepare an OpenMM simulation from a serialized thermodynamics parameter set'''
    # clear forces added from another set of thermodynamic parameters to avoid thermostat/barostat "bleedover"
    preexisting_force_indices : list[int] = sorted(
        (i for i, force in enumerate(system.getForces()) if force_added_by_polymerist(force)),
        reverse=True, # must remove greatest-to-least, as any lower-first would modify subsequent indices and cause offset
    )
    for index in preexisting_force_indices:
        LOGGER.warning(f'Removing Force "{system.getForce(index).getName()}" ({index=}) encountered from prior System setup')
        system.removeForce(index)

    # register new custom Forces (if any) from provided thermodynamic parameters
    for force in thermo_params.forces():
        LOGGER.info(f'Registering new Force "{force.getName()}" to System to enforce chosen ensemble ({thermo_params.ensemble})')
        force.setForceGroup(_POLYMERIST_FORCE_GROUP)
        system.addForce(force)
    impose_unique_force_groups(system) # NOTE: turns out to be necessary to attain energy contribution separation (can't just relabel after the simulation)

    # initialize Simulation
    simulation = Simulation(
        topology=topology,
        system=system,
        integrator=thermo_params.integrator(time_step),
    )
    
    ## set Positions
    if positions is not None:
        # TODO : add position type/shape checking
        LOGGER.info('Setting positions in Context')
        simulation.context.setPositions(positions) # set positions if provided
        
    ## set State
    state = load_state_flexible(state)
    if (state is not None):
        LOGGER.info('Setting simulation state')
        simulation.context.setState(state)

    ## ensure changes took, preserving State as necessary
    simulation.context.reinitialize(preserveState=True) # TOSELF : unclear whether this is necessary, redundant, or in fact harmful

    return simulation

def initialize_simulation_and_files(
        out_dir : Path,
        prefix : str,
        sim_params : SimulationParameters,
        topology : Topology,
        system : System,
        positions : Optional[Quantity]=None,
        state : Optional[StateLike]=None,
    ) -> tuple[Simulation, SimulationPaths]:
    '''Create simulation, bind Reporters, and update simulation Paths with newly-generated files'''
    sim_paths = SimulationPaths.from_dir_and_parameters(out_dir, prefix, sim_params, touch=True)
    simulation = simulation_from_thermo(
        topology,
        system,
        sim_params.thermo_params,
        time_step=sim_params.integ_params.time_step,
        positions=positions,
        state=state,
    )
    for reporter in sim_params.reporter_params.prepare_reporters(report_interval=sim_params.integ_params.report_interval):
        simulation.reporters.append(reporter) # add reporters to simulation instance

    return simulation, sim_paths