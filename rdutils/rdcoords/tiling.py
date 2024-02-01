'''For generating periodically-tiled topologies from RDKit Mols'''

import numpy as np

from rdkit import Chem
from rdkit.Chem import rdGeometry, rdMolTransforms

from ..rdtypes import RDMol
from ...maths import lattices
from ...maths.linearalg import affine


def tile_lattice_with_rdmol(rdmol : RDMol, unit_lattice : np.ndarray, rotate_randomly : bool=True) -> RDMol:
    '''
    Generate a tiled topology of copies of an RDKit molecule, transformed to occupy the same relative positions as points on the given lattice with unit dimensions
    if "random_rotations" is True, each occurrence of the molecule will also have a random rotation applied to it. RDMol is NOT modified at any point in the procedure
    '''
    if rdmol.GetNumConformers() < 1:
        raise ValueError('Molecule must have at least one conformer to be tiled')
    
    conformer = rdmol.GetConformer(0)
    positions = conformer.GetPositions()
    centroid  = rdMolTransforms.ComputeCentroid(conformer) # = positions.mean(axis=0) # TOSELF : mean positiona dn centroid give slightly different result, but don't affect validity of tiling
    centering = rdMolTransforms.ComputeCanonicalTransform(conformer, centroid) # translation which centers the given conformer

    dists_to_centroid = np.linalg.norm(positions - centroid, axis=1)
    r_eff = dists_to_centroid.max() # effective radius of the molecule; ensures adjacent tiled copied don't intersect
    scaled_lattice = 2 * r_eff * unit_lattice # generate scaled lattice with non-unit radius

    tiled_topology = None
    for point in scaled_lattice:
        translation = affine.xyzTrans(*point) # translation to the current point on the scaled lattice
        rotation    = affine.randRot(about_x=rotate_randomly, about_y=rotate_randomly, about_z=rotate_randomly) # rotation about either all axes or the identity matrix
        transform = translation @ rotation @ centering # NOTE : (right-to-left) order matters here!! Must first move to origin, then rotate, and then translate to latice point 

        clone = Chem.Mol(rdmol) # make copy to avoid modifying ariginal
        rdMolTransforms.TransformConformer(clone.GetConformer(0), transform) # apply transform to 
        
        if tiled_topology is None:
            tiled_topology = clone
        else:
            tiled_topology = Chem.CombineMols(tiled_topology, clone)

    return tiled_topology