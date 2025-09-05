'''For obtaining, scaling, and manipulating box vectors for Topologies'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Annotated, Literal, Union
import numpy.typing as npt

import numpy as np

from openmm.unit import Quantity
from pint import (
    Unit as PintUnit,
    Quantity as PintQuantity,
)

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
    '''
    Convert 3-vector of XYZ box dimensions into monoclinic box vectors in 3x3 diagonal reduced form
    
    Parameters
    ----------
    xyz : VectorQuantity
        3-element vector with associated units whose components
        indicate the length of the box in the x, y, and z directions

    Returns
    -------
    box_vectors : BoxVectorsQuantity
        A 3x3 matrix with associated units (same units as input vector)
        representing the monoclinic box vectors in reduced form
    '''
    assert(xyz.shape) == (3,)
    return np.diag(xyz.magnitude) * xyz.units # convert to diagonal matrix

def get_topology_bbox(offtop : Topology) -> BoxVectorsQuantity:
    '''
    Get the tight bounding box of an OpenFF Topology
    
    Parameters
    ----------
    offtop : Topology
        An OpenFF Topology

    Returns
    -------
    box_vectors : BoxVectorsQuantity
        The unit-aware box vectors representing the smallest monoclinic
        bounding box which contains all atoms in the Topology
    '''
    return xyz_to_box_vectors(offtop.get_positions().ptp(axis=0))

def _pad_box_vectors_unitless(box_vectors_mag : BoxVectors, pad_vec_mag : Vector) -> BoxVectors:
    '''Pad box vectors along each axis by a specified amount (given by each scalar component of a padding_vector)'''
    assert(pad_vec_mag.shape == (3,))

    lengths = np.linalg.norm(box_vectors_mag, axis=1) # get array of lengths of each box vector. NOTE : this NEEDS to be done on the magnitude array (i.e. can't have units) to avoid typing error
    scale_matr = np.eye(3) + 2*np.diag(pad_vec_mag / lengths) # diagonal scaling matrix, scales each vector row up to l + 2*d

    return scale_matr @ box_vectors_mag # compute matrix product to scale each row. TODO : maybe convert this to row-wise product for efficiency?
    
@allow_openmm_units
def pad_box_vectors(box_vectors : BoxVectorsQuantity, pad_vec : VectorQuantity) -> BoxVectorsQuantity:
    '''
    Pad an array of box vectors by some fixed distance from each face
    Padding amount is set for each axis separately (i.e. x, y, and z)

    Parameters
    ----------
    box_vectors : BoxVectorsQuantity
        The unit-aware box vectors array to pad
    pad_vec : VectorQuantity
        The padding vector specifying the amount to pad each pair of faces along each axis
        
        E.g. a vector of np.array([1, 2, 3])*nanometer would pad the faces perpendicular to
        the x axis by 1 nm, perpendicular to the y axis by 2 nm, and the z-axis by 3 nm
        The lengths of the edges of the box would increase by 2, 4, and 6 nm along the x, y, and z axes, respectively

    Returns
    -------
    box_vectors : BoxVectorsQuantity
        The padded box vectors
    '''
    return _pad_box_vectors_unitless(box_vectors.magnitude, pad_vec.m_as(box_vectors.units)) * box_vectors.units

def pad_box_vectors_uniform(box_vectors : BoxVectorsQuantity, pad_amount : Quantity) -> BoxVectorsQuantity:
    '''
    Pad an array of box vectors by a fixed distance from each face of the box
    
    Parameters
    ----------
    box_vectors : BoxVectorsQuantity
        The unit-aware box vectors array to pad
    pad_amount : Quantity
        The distance to pad out from each face of the box
        
        For a box with initial dimensions (Lx, Ly, Lz), the dimensions of the box returned
        after padding here will be (Lx + 2*pad_amount, Ly + 2*pad_amount, Lz + 2*pad_amount)

    Returns
    -------
    box_vectors : BoxVectorsQuantity
        The padded box vectors
    '''
    return pad_box_vectors(box_vectors, np.ones(3)*pad_amount)

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


# PROPERTIES OF BOX VECTORS
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

def encloses_box_vectors(outer_box_vectors : BoxVectorsQuantity, inner_box_vectors : BoxVectorsQuantity) -> bool:
    '''Check whether the unit cell defined by the inner box vectors is completely contained within the outer box vectors'''
    outer_vecs = box_vectors_flexible(outer_box_vectors) # DEVNOTE: done to avoid case work for missing units, wrong quantity type, etc.
    inner_vecs = box_vectors_flexible(inner_box_vectors)
    
    vec_unit : PintUnit = outer_vecs.units
    # vec_unit : PintUnit = inner_vecs.units # DEVNOTE: either choice is fine if unit dimensions are compatible; picked outer arbitrarily
    outer_vecs_unitless = outer_vecs.m_as(vec_unit)
    inner_vecs_unitless = inner_vecs.m_as(vec_unit)
    
    # Sufficient condition for containment: if the outer unit cell contains the inner, then all points in the inner unit cell,
    # including namely the inner basis vectors, are convex combinations of the outer box basis vectors
    # i.e. there exists a matrix C such that N = O @ C, where N is the inner basis, O the outer, and all elements of C are in [0, 1]
    # Equivalently, one can think of O^-1 as taking the points it contains into the standard unit cube, all of whose coordinates are in [0, 1]
    change_of_basis_matrix = np.linalg.inv(outer_vecs_unitless) @ inner_vecs_unitless
    
    return np.all((change_of_basis_matrix >= 0.0) & (change_of_basis_matrix <= 1.0)) 

