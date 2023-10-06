'''For packing solvents into Topology boxes using packmol'''

import logging
LOGGER = logging.Logger(__name__)

from typing import Optional, Union

import numpy as np
from openmm.unit import Quantity, AVOGADRO_CONSTANT_NA
from openff.toolkit import Topology, Molecule, ForceField
from openff.interchange.components import _packmol as packmol

from . import boxvectors, physprops


def pack_topology_with_solvent(offtop : Topology, solvent : Molecule, box_vecs : Union[boxvectors.VectorQuantity, boxvectors.BoxVectorsQuantity], density : Quantity, exclusion : Optional[Quantity]=None) -> Topology:
    '''Pack a Topology with a given solvent molecule to a given box size with desired density
    The dimensions of the desired box_dimensions '''
    if box_vecs.ndim == 1:
        box_vecs = boxvectors.xyz_to_box_vectors(box_vecs) # convert to 2D box vectors if only XYZ dimensions are provided

    # obtain bounding box (plus exclusion) for topology
    min_box_vecs = boxvectors.get_topology_bbox(offtop)
    if exclusion is not None:
        min_box_vecs = boxvectors.pad_box_vectors_uniform(min_box_vecs, exclusion)

    if not np.all(box_vecs >= min_box_vecs): # TODO : change this to only check XYZ dimensions
        raise boxvectors.BoxVectorError(f'Desired box dimensions ({box_vecs}) are smaller than minimum excluded Topology dimensions ({min_box_vecs})')

    # assign properties from desired box vectors if size checks pass
    box_vol = boxvectors.get_box_volume(box_vecs, units_as_openm=True)
    N = physprops.num_mols_in_box(solvent.to_rdkit(), box_vol=box_vol, density=density)

    # packing waters into topology and saving
    LOGGER.info(f'Solvating {box_vol} Topology with {N} {solvent.name} molecules to density of {density}')
    packed_top = packmol.pack_box(
        [solvent],
        [N],
        offtop,
        box_vectors=box_vecs, 
        box_shape=packmol.UNIT_CUBE,
        center_solute='BRICK'
    )
    packed_top.box_vectors = box_vecs

    return packed_top