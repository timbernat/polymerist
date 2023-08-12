'''Boilerplate for setting up OpenMM Simulations and related files'''

import logging
LOGGER = logging.getLogger(__name__)

from typing import Iterable, Optional, Union
from pathlib import Path

from openmm.app import Simulation
from openmm.app import PDBReporter, DCDReporter, StateDataReporter, CheckpointReporter

Reporter = Union[PDBReporter, DCDReporter, StateDataReporter, CheckpointReporter] # for clearer typehinting
TRAJ_REPORTERS = { # index output formats of reporters by file extension
    '.dcd' : DCDReporter,
    '.pdb' : PDBReporter
}

from .records import SimulationPaths, SimulationParameters


def prepare_simulation_paths(output_folder : Path, output_name : str, sim_params : SimulationParameters) -> SimulationPaths:
    '''Takes a Simulation object, performs energy minimization, and runs simulation for specified number of time steps
    Recording PDB frames and the specified property data to CSV at the specified frequency'''
    # creating paths to requisite output files
    prefix = f'{output_name}{"_" if output_name else ""}'
    sim_paths_out = output_folder / f'{prefix}sim_paths.json'
    sim_paths = SimulationPaths(
        sim_params=output_folder / f'{prefix}sim_parameters.json',
        trajectory=output_folder / f'{prefix}traj.{"dcd" if sim_params.binary_traj else "pdb"}',
        state_data=output_folder / f'{prefix}state_data.csv',
        checkpoint=output_folder / f'{prefix}checkpoint.{"xml" if sim_params.save_state else "chk"}',
    )
    sim_paths.to_file(sim_paths_out)
    LOGGER.info(f'Generated simulation record files at {sim_paths_out}')

    return sim_paths

def prepare_simulation_reporters(sim_paths : SimulationPaths, sim_params : SimulationParameters) ->  tuple[Reporter]:
    '''Takes a Simulation object, performs energy minimization, and runs simulation for specified number of time steps
    Recording PBD frames and the specified property data to CSV at the specified frequency'''

    # for saving pdb frames and reporting state/energy data - NOTE : all file paths must be stringified for OpenMM
    TrajReporter = TRAJ_REPORTERS[sim_paths.trajectory.suffix] # look up reporter based on the desired trajectory output file format
    
    traj_rep  = TrajReporter(file=str(sim_paths.trajectory) , reportInterval=sim_params.record_freq)  # save frames at the specified interval
    check_rep = CheckpointReporter(str(sim_paths.checkpoint), reportInterval=sim_params.record_freq, writeState=sim_params.save_state)
    state_rep = StateDataReporter(str(sim_paths.state_data) , reportInterval=sim_params.record_freq, totalSteps=sim_params.num_steps, **sim_params.reported_state_data)

    return (traj_rep, check_rep, state_rep)

def config_simulation(simulation : Simulation, reporters : Iterable[Reporter], checkpoint_path : Optional[Path]=None) -> None:
    '''Takes a Simulation object, adds data Reporters, saves an initial checkpoint, and performs energy minimization'''
    for rep in reporters:
        simulation.reporters.append(rep) # add any desired reporters to simulation for tracking

    if checkpoint_path is not None:
        simulation.saveCheckpoint(str(checkpoint_path)) # save initial minimal state to simplify reloading process
        LOGGER.info(f'Saved simulation checkpoint at {checkpoint_path}')

def run_simulation(simulation : Simulation, sim_params : SimulationParameters, output_folder : Path, output_name : str) -> None:
    '''
    Initializes an OpenMM simulation from a SMIRNOFF Interchange in the desired ensemble
    Creates relevant simulation files, generates Reporters for state, checkpoint, and trajectory data,
    performs energy minimization, then integrates the trajectory for the desired number of steps
    '''
    sim_paths = prepare_simulation_paths(output_folder, output_name, sim_params)
    reporters = prepare_simulation_reporters(sim_paths, sim_params)
    sim_params.to_file(sim_paths.sim_params) # TOSELF : this is not a parameters checkpoint file UPDATE, but rather the initial CREATION of the checkpoint file
    config_simulation(simulation, reporters, checkpoint_path=sim_paths.checkpoint)

    LOGGER.info('Performing energy minimization')
    simulation.minimizeEnergy()

    LOGGER.info(f'Integrating {sim_params.total_time} OpenMM sim at {sim_params.temperature} and {sim_params.pressure} for {sim_params.num_steps} steps')
    simulation.step(sim_params.num_steps)
    LOGGER.info('Simulation integration completed successfully')
