'''For obtaining, scaling, and manipulating box vectors for Topologies'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Annotated, Literal, TypeAlias, Union
import numpy.typing as npt

import numpy as np

from openmm.unit import Quantity
from pint import Quantity as PintQuantity

from openff.toolkit import Topology
from openff.interchange.components._packmol import _box_vectors_are_in_reduced_form

from .unitsys import allow_openmm_units, openff_to_openmm


# CUSTOM TYPES FOR CLARITY, ESPECIALLY WITH UNITS
Vector = Annotated[npt.NDArray[np.generic], Literal[3]] # a 3x1 vector
VectorQuantity = Union[Quantity, Vector] # 3x1 vector with associated units

BoxVectors = Annotated[npt.NDArray[np.generic], Literal[3, 3]] # a 3x3 box-vector matrix (in reduced form), with each row representing a single unit cell vector
BoxVectorsQuantity = Union[Quantity, BoxVectors] # 3x3 box vectors with associated units

class BoxVectorError(Exception):
    '''Raised when a provided set of box vectors is invalid (for whatever reason)'''
    pass


# OBTAINING AND SCALING BOX VECTORS
@allow_openmm_units
def xyz_to_box_vectors(xyz : VectorQuantity) -> BoxVectorsQuantity:
    '''Convert 3-vector of XYZ box dimensions into monoclinic box vectors in 3x3 diagonal reduced form'''
    assert(xyz.shape) == (3,)
    return np.diag(xyz.magnitude) * xyz.units # convert to diagonal matrix

def get_topology_bbox(offtop : Topology) -> BoxVectorsQuantity:
    '''Get the tight bounding-box XYZ dimensions of a Topology'''
    xyz = offtop.get_positions().ptp(axis=0)
    return xyz_to_box_vectors(xyz)

def _pad_box_vectors_unitless(box_vectors_mag : BoxVectors, pad_vec_mag : Vector) -> BoxVectors:
    '''Pad box vectors along each axis by a specified amount (given by each scalar component of a padding_vector)'''
    assert(pad_vec_mag.shape == (3,))

    lengths = np.linalg.norm(box_vectors_mag, axis=1) # get array of lengths of each box vector. NOTE : this NEEDS to be done on the magnitude array (i.e. can't have units) to avoid typing error
    scale_matr = np.eye(3) + 2*np.diag(pad_vec_mag / lengths) # diagonal scaling matrix, scales each vector row up to l + 2*d

    return scale_matr @ box_vectors_mag # compute matrix product to scale each row. TODO : maybe convert this to row-wise product for efficiency?
    
@allow_openmm_units
def pad_box_vectors(box_vectors : BoxVectorsQuantity, pad_vec : VectorQuantity) -> BoxVectorsQuantity:
    '''Pad each box vector on either side by a fixed distance given by a component of a padding vector'''
    return _pad_box_vectors_unitless(box_vectors.magnitude, pad_vec.m_as(box_vectors.units)) * box_vectors.units

def pad_box_vectors_uniform(box_vectors : BoxVectorsQuantity, pad_amount : Quantity) -> BoxVectorsQuantity:
    '''Padd all box vectors by the same fixed distance'''
    return pad_box_vectors(box_vectors, np.ones(3)*pad_amount)


# OBTAINING VOLUMES
def _get_box_volume_unitless(box_vectors_mag : BoxVectors) -> float:
    '''Compute volume of a box given a set of 3x3 box vectors'''
    return np.abs(np.linalg.det(box_vectors_mag)) # return unsigned determinant

def get_box_volume(box_vectors : BoxVectorsQuantity, units_as_openm : bool=True) -> Quantity:
    '''Compute volume of a box given a set of 3x3 box vectors'''
    assert(_box_vectors_are_in_reduced_form(box_vectors))
    box_vol = _get_box_volume_unitless(box_vectors.magnitude) * box_vectors.units**3

    if units_as_openm:
        return openff_to_openmm(box_vol)
    return box_vol


# CONVERSIONS
@allow_openmm_units
def box_vectors_flexible(box_vecs : Union[VectorQuantity, BoxVectorsQuantity]) -> BoxVectorsQuantity:
    '''Allows for passing XYZ box dimensions '''
    if not isinstance(box_vecs, PintQuantity):
        raise TypeError('Box vectors passed have no associated units')
    
    if not isinstance(box_vecs.magnitude, np.ndarray):
        raise TypeError('Must pass array-like Quantity as box vector')
    
    if box_vecs.ndim == 1:
        return xyz_to_box_vectors(box_vecs)
    return box_vecs
