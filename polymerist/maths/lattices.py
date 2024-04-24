'''Utilities for generating periodic unit lattices'''

from typing import Callable, ClassVar, Iterable, Optional, Type
from dataclasses import dataclass, field

import numpy as np
from itertools import product as cartesian_product

from ..genutils.typetools.numpytypes import Shape, N, M, DType


# COORDINATE ARRAY OPERATIONS
def mean_coord(coords : np.ndarray[Shape[M, N], DType]) -> np.ndarray[Shape[N], DType]:
    '''The average (center-of-mass) coordinate of a vector of coordinates'''
    return coords.mean(axis=0)
center_of_mass = COM = mean_coord

def dists_to_point(coords : np.ndarray[Shape[M, N], DType], point : np.ndarray[Shape[N], DType], norm_order : Optional[int]=None) -> np.ndarray[Shape[M], DType]:
    '''The distance between each point in a coordinate array and a single arbitrary point'''
    return np.linalg.norm(coords - point, ord=norm_order, axis=1)

def dists_to_centroid(coords : np.ndarray[Shape[M, N], DType], norm_order : Optional[int]=None) -> np.ndarray[Shape[M], DType]:
    '''The distance of each coordinate in an array of coordinates to the coordinates' centroid'''
    return dists_to_point(coords, point=mean_coord(coords), norm_order=norm_order)

def bounding_box_extrema(coords : np.ndarray[Shape[M, N], DType]) -> np.ndarray[Shape[2, N], DType]:
    '''The minimal and maximal coordinates of the tight bounding box around an array of coordinate'''
    return coords.min(axis=0), coords.max(axis=0)

def bounding_box_dims(coords : np.ndarray[Shape[M, N], DType]) -> np.ndarray[Shape[N], DType]:
    '''The tight bounding box side lengths around an array of coordinates'''
    return coords.ptp(axis=0)


# INTEGER LATTICES (UNIT BASIS SUPERCELLS)
def generate_int_lattice(*dims : Iterable[int]) -> np.ndarray[Shape[N, 3], int]:
    '''Generate all N-D coordinates of points on a integer lattice with the sizes of all D dimensions given'''
    return np.fromiter(
        iter=cartesian_product(*[
            range(d)
                for d in dims
        ]),
        dtype=np.dtype((int, len(dims)))
    )


# LATTICE VECTORS AND PARAMETERS
@dataclass
class LatticeParameters: # TODO : add support for OpenMM units, add compat w/ lammpstools.lammpseval
    '''For storing the lengths of and the angles between the 3 lattice vectors of a crystallographic unit cell'''
    a : float
    b : float
    c : float

    alpha : float = field(default=np.pi / 2) # make cell orthorhombic by default
    beta  : float = field(default=np.pi / 2) # make cell orthorhombic by default
    gamma : float = field(default=np.pi / 2) # make cell orthorhombic by default

    _LENGTH_ATTRS : ClassVar[tuple[str]] = ('a', 'b', 'c')
    _ANGLE_ATTRS  : ClassVar[tuple[str]] = ('alpha', 'beta', 'gamma')

    # VALIDATION
    def __post_init__(self) -> None:
        '''Validating initialized values'''
        lens_positive, bad_len_name = self.lengths_are_positive
        if not lens_positive:
            raise ValueError(f'{self.__class__.__name__} unit cell lengths must be positive (passed invalid value {bad_len_name}={getattr(self, bad_len_name)})')
        self.reduce_angles()

    @property
    def lengths_are_positive(self) -> tuple[bool, Optional[str]]:
        '''Check whether all axial lengths are well-defined (i.e. positive)'''
        for len_attr in self._LENGTH_ATTRS:
            if (getattr(self, len_attr) <= 0):
                return False, len_attr
        else:
            return True, None
        
    def reduce_angles(self) -> None:
        '''Ensure all angles are within the principal [0, 2*pi) interval'''
        for angle, angle_attr in zip(self.angles, self._ANGLE_ATTRS):
            setattr(self, angle_attr, angle % (2*np.pi)) # ensure angles are within [0, 2*pi)

    # PROPERTY VIEWS
    @property
    def axial_lengths(self) -> np.ndarray[Shape[3], float]:
        '''View of just the axial lengths'''
        return np.array([getattr(self, len_attr) for len_attr in self._LENGTH_ATTRS])
    lengths = axial_lengths

    @property
    def angles(self) -> np.ndarray[Shape[3], float]:
        '''Property alias of angles() method for convenience'''
        return np.array([getattr(self, len_attr) for len_attr in self._ANGLE_ATTRS])
    
    def axial_angles(self, in_degrees : bool=False) -> np.ndarray[Shape[3], float]:
        '''View of just the angles'''
        if in_degrees:
            return np.rad2deg(self.angles)
        return self.angles
    
    @property
    def volume(self) -> float:
        '''The volume of the unit cell (in arbitrary units)'''
        return np.linalg.det(self.lattice_vectors)

    # LATTICE VECTOR METHODS
    @classmethod
    def from_lattice_vectors(cls, vectors : np.ndarray[Shape[3,3], float]) -> 'LatticeParameters':
        '''Obtain axial lengths and inter-vector angles from a matrix of lattice vectors'''
        axial_lengths = np.linalg.norm(vectors, axis=1)
        normed_vectors = vectors / axial_lengths[:, None] # normalize along rows

        BCA = np.roll(normed_vectors, -1, axis=0) # lattice vectors cycled backward
        CAB = np.roll(normed_vectors,  1, axis=0) # lattice vectors cycled forward
        cycled_dots = np.sum(BCA * CAB, axis=1) # row-wise dot product between cycled rows, giving normalized [B.C, C.A, and A.B] as result
        angles = np.arccos(cycled_dots) # invert cosine from dot product expression to recover angles

        return cls(*axial_lengths, *angles)
    
    @staticmethod
    def _lattice_vector_classmethod_factory(lattice_vectors : np.ndarray[Shape[3,3]], lattice_name : str, prefix : str='create') -> Callable[[Type['LatticeParameters'], float], 'LatticeParameters']:
        '''Factory function for generating LatticeParameter initializers based on pre-defined lattice vectors'''
        def generic_lattice_generator(cls, dim : float=1.0) -> LatticeParameters:
            return cls.from_lattice_vectors(dim * lattice_vectors)
        
        # modify signature of generic to mask its dynamic creation
        fn_name = f'{prefix}{"_" if prefix else ""}{lattice_name.lower()}'
        generic_lattice_generator.__name__ = fn_name
        generic_lattice_generator.__doc__  = f'Generate lattice parameters for a {lattice_name.lower()} box (with optional uniform scaling)'
        generic_lattice_generator.__qualname__ = f'{LatticeParameters.__qualname__}.{fn_name}'

        return generic_lattice_generator

    def to_lattice_vectors(self) -> np.ndarray[Shape[3,3], float]:
        '''
        The restricted lattice vectors corresponding to the lattice parameters,
        where vector A lies along the x-axis and vector B is in the xy-plane
        
        Vectors are returned as a 3x3 matrix, where each row represents a lattice vector
        '''
        # formulas adapted from https://www.aflowlib.org/prototype-encyclopedia/triclinic_lattice.html 
        ax = self.a * 1.0                # vector A lies purely along the x-axis
        bx = self.b * np.cos(self.gamma) #    projection onto the x-axis unit vector (i.e. the normalized A vector) via the dot product
        by = self.b * np.sin(self.gamma) # perpendicular onto the x-axis unit vector (i.e. the normalized A vector) via the cross product
        cx = self.c * np.cos(self.beta)  #    projection onto the x-axis unit vector (i.e. the normalized A vector) via the dot product
        cy = self.c * (np.cos(self.alpha) - np.cos(self.beta)*np.cos(self.gamma)) / np.sin(self.gamma) # projection of C onto the y-axis unit vector (via the vector triple product C*(AxB)xA)
        cz = np.sqrt(self.c**2 - cx**2 - cy**2) # remaining z-axis contribution to C must all be along z-axis (algebraically simpler this way)

        return np.array([
            [ax , 0.0, 0.0],
            [bx , by , 0.0],
            [cx , cy , cz ],
        ])

    @property
    def lattice_vectors(self) ->np.ndarray[Shape[3,3], float]:
        '''Property alias of to_lattice_vectors() method for convenience'''
        return self.to_lattice_vectors()
    
## COMMON UNIT LATTICE VECTORS FOR REFERENCE (WITH DYNAMIC REGISTRATION)
COMMON_LATTICE_VECTORS : dict[str, np.ndarray[Shape[3, 3], float]] = {
    'UNIT_CUBE' : np.asarray([
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ]),
    'RHOMBIC_DODECAHEDRON_XY_SQR' : np.array([
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1/2, 1/2, np.sqrt(2.0) / 2.0],
    ]),
    'RHOMBIC_DODECAHEDRON_XY_HEX' : np.array([
        [1.0, 0.0, 0.0],
        [1/2, np.sqrt(3.0) / 2.0, 0.0],
        [1/2, np.sqrt(3.0) / 6.0, np.sqrt(6.0) / 3.0],
    ]),
    'TRUNCATED_OCTAHEDRON' : np.array([
        [ 1.0, 0.0, 0.0],
        [ 1/3, 2.0 * np.sqrt(2.0) / 3.0, 0.0],
        [-1/3, np.sqrt(2.0) / 3.0, np.sqrt(6.0) / 3.0],
    ]),
} # adapted from https://manual.gromacs.org/current/reference-manual/algorithms/periodic-boundary-conditions.html

for lattice_name, lattice_vectors in COMMON_LATTICE_VECTORS.items():
    globals()[lattice_name] = lattice_vectors # dynamically register common lattice vectors at module level
    
    lattice_generator = LatticeParameters._lattice_vector_classmethod_factory(lattice_vectors, lattice_name=lattice_name)
    setattr(LatticeParameters, lattice_generator.__name__, classmethod(lattice_generator))
