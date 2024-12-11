'''For interfacing between OpenFF and OpenMM representations of Topologies and other MD primitives'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Optional, Union

from pathlib import Path

from openff.toolkit import ForceField
from openff.interchange import Interchange
from openff.toolkit.topology import Topology as OFFTopology

from openmm import System
from openmm.app import Topology as OMMTopology
from openmm.unit import Quantity

from ..unitsys import openff_to_openmm
from .. import FFDIR
from ..boxvectors import box_vectors_flexible, VectorQuantity, BoxVectorsQuantity


def forcefield_flexible(forcefield : Union[ForceField, str, Path]) -> ForceField: # DEV: consider deprecating
    '''For making forcefield input to other functions more flexible (can accept a literal ForceField, a string name, or a Path to the forcefield)'''
    if isinstance(forcefield, ForceField): 
        return forcefield
    
    if isinstance(forcefield, str): # NOTE : order here matters, contiguous with subsequent Path clause (DO NOT put this after it)
        forcefield = Path(forcefield)

    if isinstance(forcefield, Path):
        if FFDIR in forcefield.parents:
            ff_path = forcefield #  if the Path already contains the forcefield directory, presume it to be a full forcefield path
        else:
            ff_path = FFDIR / forcefield # otherwise, presume it to be the name of a forcefield

        if ff_path.suffix != '.offxml':
            ff_path = ff_path.with_name(f'{ff_path.name}.offxml') # TODO : find more general way to handle the presence of dots in FF names (e.g. openff-2.0.0) in conjunction with suffixes

        return ForceField(ff_path)

def openff_topology_to_openmm(
            offtop : OFFTopology,
            forcefield : Union[ForceField, str, Path],
            box_vecs : Optional[Union[VectorQuantity, BoxVectorsQuantity]]=None,
            combine_nonbonded_forces : bool=False,
            add_constrained_forces : bool=False
        ) -> tuple[OMMTopology, System, Quantity]:
    '''Converts an OpenFF Topology to an OpenMM Topology, System, and Positions'''
    if box_vecs is not None:
        offtop.box_vectors = box_vectors_flexible(box_vecs)

    forcefield = forcefield_flexible(forcefield) 
    ic = Interchange.from_smirnoff(forcefield, offtop, charge_from_molecules=[mol for mol in set(offtop.molecules)]) # WARNING : will call toolkits wrappers for charge assignment if Topology is not fully charged!

    ommtop = ic.to_openmm_topology()
    ommsys = ic.to_openmm(combine_nonbonded_forces=combine_nonbonded_forces, add_constrained_forces=add_constrained_forces)
    ommpos = openff_to_openmm(ic.positions)

    return ommtop, ommsys, ommpos