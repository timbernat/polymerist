'''For packing solvents into Topology boxes using packmol'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import logging
LOGGER = logging.getLogger(__name__)
from typing import Union

from openmm.unit import Quantity
from openff.toolkit import Topology, Molecule
from openff.interchange.components import _packmol as packmol

from ..boxvectors import (
    VectorQuantity,
    BoxVectorsQuantity,
    BoxVectorError,
    box_vectors_flexible,
    get_topology_bbox,
    get_box_volume,
    encloses_box_vectors,
)   
from ..physprops import num_mols_in_box
from ..topology import get_largest_offmol


def pack_topology_with_solvent(
        offtop : Topology,
        solvent : Molecule,
        box_vecs : Union[VectorQuantity, BoxVectorsQuantity],
        density : Quantity,
        **kwargs,
    ) -> Topology:
    '''Pack a Topology with a given solvent molecule to a given box size with desired density
    The dimensions of the desired box_dimensions '''
    # obtain bounding box (plus exclusion) for topology
    box_vecs = box_vectors_flexible(box_vecs) # up-convert to full box vectors if only dimensions are given
    min_box_vecs = get_topology_bbox(offtop)
    if not encloses_box_vectors(box_vecs, min_box_vecs):
        raise BoxVectorError(f'Desired box vectors ({box_vecs}) are smaller than minimum possible vectors which enclose the Topology ({min_box_vecs})')

    # assign properties from desired box vectors if size checks pass
    box_vol = get_box_volume(box_vecs, units_as_openm=True)
    N = num_mols_in_box(solvent, box_vol=box_vol, density=density)

    # packing waters into topology and saving
    LOGGER.info(f'Solvating {box_vol} Topology with {N} {solvent.name} molecules to density of {density}')
    packed_top = packmol.pack_box(
        molecules=[solvent],
        number_of_copies=[N],
        solute=offtop,
        box_vectors=box_vecs, 
        box_shape=packmol.UNIT_CUBE,
        center_solute='BRICK',
        **kwargs
    )
    LOGGER.info('Packmol packing converged')

    packed_top.box_vectors = box_vecs
    LOGGER.info(f'Set solvated Topology box vectors to {box_vecs}')
    offmol = get_largest_offmol(packed_top)
    offmol.properties['solvent'] = solvent.name

    return packed_top