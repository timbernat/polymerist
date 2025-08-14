'''For running OpenMM simulations and extracting information from them'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import logging
LOGGER = logging.getLogger(__name__)

from typing import Optional, Union
from pathlib import Path
from numpy import ndarray

from openmm import System
from openmm.app import Simulation, Topology
from openmm.unit import Quantity

from .parameters import SimulationParameters
from .preparation import initialize_simulation_and_files, StateLike
from .serialization.system import serialize_system
from .serialization.topology import serialize_topology_from_simulation
from .serialization.paths import SimulationPaths


def run_simulation_schedule(
        working_dir : Path,
        schedule : dict[str, SimulationParameters],
        init_top : Topology,
        init_sys : System,
        init_pos : ndarray,
        init_state : Optional[StateLike]=None,
        return_history : bool=False,
    ) -> Optional[dict[str, tuple[Simulation, SimulationPaths]]]:
    '''Run several OpenMM simulations in series, based on an initial set of OpenMM objects and a "schedule" consisting of a sequence of named parameter sets'''
    if not isinstance(init_pos, Quantity):
        raise TypeError('Positions must have associated OpenMM units') # TODO : provide more robust check for this
    
    working_dir.mkdir(exist_ok=True)
    history = {}

    num_steps = len(schedule)
    for i, (step_name, sim_params) in enumerate(schedule.items()): # TOSELF : may want to shift to an explicitly-ordered map (e.g. collections OrderedDict); Python 3.6+ preserves order, but wouldn't hurt to be extra safe
        if i == 0:
            ommtop = init_top
            ommsys = init_sys
            ommpos = init_pos
            ommstate = init_state
        else:
            ommtop = simulation.topology
            ommsys = simulation.system
            ommpos = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)
            ommstate = None # DEVNOTE: specifically DON'T want to bleed positions, velocities, etc. between simulations with potentially-different thermodynamic parameters

        LOGGER.info(f'Initializing simulation {i + 1}/{num_steps} ("{step_name}")')
        simulation, sim_paths = initialize_simulation_and_files(
            out_dir=working_dir / step_name,
            prefix=step_name,
            sim_params=sim_params,
            topology=ommtop,
            system=ommsys,
            positions=ommpos,
            state=ommstate,
        )
        history[step_name] = {
            'simulation' : simulation,
            'paths'      : sim_paths,
        }
        
        LOGGER.info(f'Performing energy minimization (initial PE = {simulation.context.getState(getEnergy=True).getPotentialEnergy()})')
        simulation.minimizeEnergy()
        LOGGER.info(f'Energy successfully minimized (final PE = {simulation.context.getState(getEnergy=True).getPotentialEnergy()})')

        serialize_topology_from_simulation(sim_paths.topology_path, simulation) # TODO : make this consistent with rest of Path output
        LOGGER.info(f'Saved energy-minimized Simulation Topology at {sim_paths.topology_path}')

        serialize_system(sim_paths.system_path, simulation.system)
        LOGGER.info(f'Saved serialized Simulation System at {sim_paths.system_path}')

        LOGGER.info(f'Integrating {sim_params.integ_params.total_time} OpenMM Simulation for {sim_params.integ_params.num_steps} steps')
        simulation.step(sim_params.integ_params.num_steps)
        LOGGER.info('Simulation integration completed successfully')
        LOGGER.info('') # add whitespace between simulations to act as a "palate cleanser"

    if return_history:
        return history