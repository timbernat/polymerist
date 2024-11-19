'''Boilerplate for setting up OpenMM Simulations and related files'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import logging
LOGGER = logging.getLogger(__name__)

from typing import Optional, Union
from pathlib import Path

from openmm import State, System, XmlSerializer
from openmm.app import Simulation, Topology
from openmm.unit import Quantity

# from .records import SimulationParameters
from .serialization import SimulationPaths
from .thermo import EnsembleFactory
from .parameters import ThermoParameters, SimulationParameters
from .forcegroups import impose_unique_force_groups


def load_state_flexible(state : Optional[Union[str, Path, State]]=None) -> Optional[State]:
    '''Allows one to flexibly load an OpenMM state, either from a State object or file-like object'''
    if isinstance(state, State) or (state is None):
        state = state
    else:
        if isinstance(state, Path):
            state_path = state
        elif isinstance(state, str):
            state_path = Path(state)
        # TODO : add support for load from opened file
        else:
            raise TypeError('State can only be loaded from pathlike object') 
        
        try:
            with state_path.open('r') as state_file:
                state = XmlSerializer.deserialize(state_file.read())
        except ValueError:
            state = None
    
    if state is None:
        LOGGER.warning('No valid State/State file provided, initializing State as None')
    return state

def simulation_from_thermo(topology : Topology, system : System, thermo_params : ThermoParameters, time_step : Quantity, state : Optional[Union[str, Path, State]]=None) -> Simulation:
    '''Prepare an OpenMM simulation from a serialized thermodynamics parameter set'''
    ens_fac = EnsembleFactory.from_thermo_params(thermo_params)
    if (extra_forces := ens_fac.forces()): # check if any extra forces are present
        for force in extra_forces:
            system.addForce(force) # add forces to System BEFORE creating Simulation to avoid having to reinitialize the Conext to preserve changes 
            LOGGER.info(f'Added {force.getName()} Force to System')
    impose_unique_force_groups(system)

    simulation = Simulation(
        topology=topology,
        system=system,
        integrator=ens_fac.integrator(time_step),
    )
    state = load_state_flexible(state)
    if state is not None:
        LOGGER.info('Setting simulation state')
        simulation.context.setState(state)
        simulation.context.reinitialize(preserveState=True) # TOSELF : unclear whether this is necessary, redundant, or in fact harmful

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
        LOGGER.info('Setting positions in Context')
        simulation.context.setPositions(positions) # set positions if provided

    return simulation, sim_paths