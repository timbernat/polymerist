'''For generating periodically-tiled topologies from RDKit Mols'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import numpy as np

from rdkit import Chem
from rdkit.Chem import rdGeometry, rdMolTransforms, Mol

from ...genutils.typetools.numpytypes import Shape, N
from ...maths.linearalg import affine
from ...maths.lattices.coordinates import Coordinates


def rdmol_effective_radius(rdmol : Chem.Mol, conf_id : int=0) -> float:
    '''Determine an effective radius of influence of a molecule such that all atoms are within this radius'''
    conformer = rdmol.GetConformer(conf_id)
    positions = Coordinates(conformer.GetPositions())
    dists_to_centroid = positions.dists_to_centroid()

    # polymerist-agnostic implementation
    # centroid  = rdMolTransforms.ComputeCentroid(conformer, ignoreHs=True) # = positions.mean(axis=0) # TOSELF : mean positiona dn centroid give slightly different result, but don't affect validity of tiling
    # dists_to_centroid = np.linalg.norm(positions - centroid, axis=1)

    return dists_to_centroid.max()

def tile_lattice_with_rdmol(rdmol : Mol, lattice_points : np.ndarray[Shape[N, 3], float], rotate_randomly : bool=True, conf_id : int=0) -> Mol:
    '''
    Generate a tiled topology of copies of an RDKit molecule, transformed to occupy the same relative positions as points on the given lattice with unit dimensions
    if "random_rotations" is True, each occurrence of the molecule will also have a random rotation applied to it. Mol is NOT modified at any point in the procedure
    '''
    if rdmol.GetNumConformers() < 1:
        raise ValueError('Molecule must have at least one conformer to be tiled')
    
    conformer = rdmol.GetConformer(conf_id)
    centering = rdMolTransforms.ComputeCanonicalTransform(conformer, ignoreHs=True) # translation which centers the given conformer

    tiled_topology = None
    for point in lattice_points:
        translation = affine.xyzTrans(*point) # translation to the current point on the lattice
        rotation    = affine.randRot(about_x=rotate_randomly, about_y=rotate_randomly, about_z=rotate_randomly) # rotation about either all axes or the identity matrix
        transform = translation @ rotation @ centering # NOTE : (right-to-left) order matters here!! Must first move to origin, then rotate, and then translate to latice point 

        clone = Chem.Mol(rdmol) # make copy to avoid modifying ariginal
        rdMolTransforms.TransformConformer(clone.GetConformer(0), transform) # apply transform to 
        
        if tiled_topology is None:
            tiled_topology = clone
        else:
            tiled_topology = Chem.CombineMols(tiled_topology, clone)

    return tiled_topology