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


# TODO : add optional_in_place??
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
        simulation.context.setPositions(positions) # set positions if provided

    return simulation, sim_paths
