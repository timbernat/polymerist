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

from .thermo import EnsembleFactory
from .parameters import ThermoParameters, SimulationParameters
from .serialization.paths import SimulationPaths
from .serialization.state import StateLike, load_state_flexible


def simulation_from_thermo(
        topology : Topology,
        system : System,
        thermo_params : ThermoParameters,
        time_step : Quantity,
        state : Optional[StateLike]=None,
    ) -> Simulation:
    '''Prepare an OpenMM simulation from a serialized thermodynamics parameter set'''
    ens_fac = EnsembleFactory.from_thermo_params(thermo_params)
    if (extra_forces := ens_fac.forces()): # check if any extra forces are present
        for force in extra_forces:
            # TODO: label injected forces with dedicated forceGroup for tracking downstream
            system.addForce(force) # add forces to System BEFORE creating Simulation to avoid having to reinitialize the Conext to preserve changes 
            LOGGER.info(f'Added {force.getName()} Force to System')

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
    if state is None:
        LOGGER.info(f'No explicit initial State supplied, defaulting to State cached in "{sim_paths.state_path}"')
        state = sim_paths.state_path

    # create simulation and add reporters
    simulation = simulation_from_thermo(
        topology,
        system,
        sim_params.thermo_params,
        time_step=sim_params.integ_params.time_step,
        state=state,
    )
    for reporter in sim_params.reporter_params.prepare_reporters(report_interval=sim_params.integ_params.report_interval):
        simulation.reporters.append(reporter) # add reporters to simulation instance

    if positions is not None:
        # TODO : add position type/shape checking
        LOGGER.info('Setting positions in Context')
        simulation.context.setPositions(positions) # set positions if provided

    return simulation, sim_paths