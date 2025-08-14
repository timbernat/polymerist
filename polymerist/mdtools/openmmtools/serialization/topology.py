'''For handling serialization of molecular topologies and positions in OpenMM format'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import logging
LOGGER = logging.getLogger(__name__)

from typing import Any, Optional, Union
from numpy.typing import NDArray

from pathlib import Path
from openmm import Vec3
from openmm.app import (
    Simulation,
    PDBFile,
    Topology as OpenMMTopology,
)
from ....genutils.fileutils.pathutils import allow_string_paths
from ....molfiles.pdb import SerialAtomLabeller


@allow_string_paths
def serialize_openmm_pdb(
        pdb_path : Path,
        topology : OpenMMTopology,
        positions : Union[NDArray, list[Vec3]],
        keep_chain_and_res_ids : bool=True,
        atom_labeller : Optional[SerialAtomLabeller]=None,
        resname_map : Optional[dict[str, str]]=None,
    ) -> None:
    '''Configure and write an Protein DataBank File from an OpenMM Topology and array of positions
    Provides options to configure atom ID numbering, residue numbering, and residue naming'''
    if atom_labeller is None:
        atom_labeller = SerialAtomLabeller()
    
    if resname_map is None:
        resname_map = {} # avoids mutable default

    # chain config
    for chain in topology.chains():
        chain.id = str(chain.id)

    # residue config
    for residue in topology.residues():
        residue.id = str(residue.id) # avoids TypeError when specifying keepIds during PDB write
        repl_res_name = resname_map.get(residue.name, None) # lookup current residue name to see if a replacement is called for
        if repl_res_name is not None:
            residue.name = repl_res_name

    # individual atom config
    if atom_labeller: # implicitly, preserves extant atom names if a labeller is not given
        for atom in topology.atoms():
            if atom.element is not None:
                atom_label = atom.element.symbol
            elif atom.name == 'EP': # "Extra Particle", as defined by SMIRNOFF spec https://openforcefield.github.io/standards/standards/smirnoff/#virtualsites-virtual-sites-for-off-atom-charges
                atom_label = atom.name # 'EP'
            else:
                continue
            
            atom.name = atom_labeller.get_atom_label(atom_label)

    # file write
    with pdb_path.open('w') as file:
        PDBFile.writeFile(topology, positions, file, keepIds=keep_chain_and_res_ids)

@allow_string_paths
def serialize_topology_from_simulation(pdb_path : Path, sim : Simulation, keep_ids : bool=False) -> None:
    '''Saves a PDB of the current state of a simulation's Topology'''
    serialize_openmm_pdb(
        pdb_path,
        topology=sim.topology,
        positions=sim.context.getState(getPositions=True).getPositions(),
        keep_chain_and_res_ids=keep_ids
    )