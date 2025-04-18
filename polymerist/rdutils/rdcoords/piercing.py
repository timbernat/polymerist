'''Methods for detecting ring piercing in molecular conformers'''

from typing import Sequence, TypeVar
N = TypeVar('N', bound=int)

import numpy as np
from scipy.spatial import Delaunay

from rdkit import Chem


def PINPRICS(
    positions : np.ndarray[tuple[N, 3], float],
    bonded_atom_idxs : Sequence[tuple[int, int]],
    ring_atom_idxs : list[int], # this one really does have to be a list (not a tuple!)
) -> tuple[tuple[int, int]]:
    '''Implements the PINPRICS algorithm (Planar Intersection Nomography for Pierced Ring Identification via Cut Sets)
    Assumes rings is oblate, i.e. faces along the minor axis and is "wider than it is tall"
    
    Parameters
    ----------
    positions : Array[[N, 3], float]
        An array of the 3D coordinates of all atoms in a molecule
    bonded_atom_idxs : Sequence[tuple[int, int]]
        A sequence of pairs of indices of all bonded atoms in the molecule
        This represents the topology on the molecules graph
    ring_atom_idxs : list[int]
        List of indices of the atoms which make up the ring to be tested for piercing
        Requirement of specifically a list (i.e. not just a Sequence etc.) is to allow for "smart" indexing in numpy

    Returns
    -------
    piercing_idxs : tuple[tuple[int, int]]
        A tuple of all pairs of atom indices which were perceived to pierce the ring
    '''
    # preliminary type recasts
    bonded_atom_idxs = np.array(bonded_atom_idxs) # convert bond indices to array to allow fancy numpy indexing - TODO: add type checks?
    if isinstance(ring_atom_idxs, tuple):
        ring_atom_idxs = list(ring_atom_idxs) # tuple screws up numpy indexing, need this to be a list (or an array)
    
    # local origin and basis for ring
    ## locate and subtract off center
    ring_positions = positions[ring_atom_idxs]
    ring_center = np.mean(ring_positions, axis=0)
    positions_centered = positions - ring_center
    ring_positions = positions_centered[ring_atom_idxs]
    
    ## identify ring normal, defined as the minimum-variance direction (assuming oblate ring)
    U, S, Vh = np.linalg.svd(ring_positions, full_matrices=False) # NOTE: SVD here requires centered ring positions (i.e. subtracting mean)
    ring_axes = Vh.T                # transpose places eigenvectors into column-order; basis guaranteed to be orthogonal, since covariance matrix is real and symmetric
    ring_normal = ring_axes[:, -1]  # NOTE: no need to sort by eigenvalue, linalg.svd already sorts in descending order by singular value

    # determine orientation of all atoms relative to ring plane:
    # *  1 <=> in front of ring
    # *  0 <=> parallel to/part of ring
    # * -1 <=> behind ring
    ring_affiliations = np.sign(np.dot(positions_centered, ring_normal)).astype(int) # NOTE: use of uncentered postiions is NOT a mistake (subtracting center alters dot product)
    ring_affiliations[ring_atom_idxs] = 0 # explicitly clamp ring atoms as being parallel (even if not strictly geometrically in-plane)
    spanning_bond_idxs = (np.prod(ring_affiliations[bonded_atom_idxs], axis=1) == -1) # given A = {-1, 0, 1}, the only pairs from AxA which multiply to -1 are (1, -1) and (-1, 1)
    spanning_pair_idxs = bonded_atom_idxs[spanning_bond_idxs, :]
    
    # interpolate spanning bonds and solve for ring plane crossings
    spanning_positions = positions_centered[spanning_pair_idxs]
    bond_starts, bond_ends = spanning_positions[:, 0, :], spanning_positions[:, 1, :]
    bond_vectors = bond_ends - bond_starts
    
    t_intersects = -np.dot(bond_starts, ring_normal) / np.dot(bond_vectors, ring_normal) # NOTE: negative sign comes from np.dot(center - start, normal), since center is now at the origin
    bond_plane_intersections = bond_starts + bond_vectors*t_intersects[:, None] # plug parameter back in (need to broadcast!) to solve for bond-plane intersections
    t_is_supported = (t_intersects >= 0) & (t_intersects <= 1) # since physically a bond only exists between its atoms, the LERP is only supported on t in [0, 1]
    
    # apply orthogonal transformation to align ring normal direction with the z-axis (the "conductor")
    ## determine projections of ring points onto ring plane
    O = planar_projector = np.eye(3, dtype=float) - np.outer(ring_normal, ring_normal)/np.inner(ring_normal, ring_normal) # !NOTE!: the "lack" of a 2 coefficient is deliberate, as this is NOT a Householder matrix (project onto -but not through- plane)
    ring_positions_planar = ring_positions @ O.T

    ## Align ring normal with +z-axis
    diff = ring_normal - np.array([0, 0, 1], dtype=float) # technically a parameter, but for polygon intersection tests we need this to be one of the standard basis vectors (picked z arbitrarily)
    H = reflector = np.eye(3, dtype=float) - 2*np.outer(diff, diff)/np.inner(diff, diff) # compute Householder matrix which reflects ring normal onto conductor (z axis unit vector here)
        
    ring_positions_planar @= H.T # apply alignment reflection
    assert np.allclose(ring_positions_planar[:, -1], 0.0) # double-check that the above manipulations have aligned the array into the XY-plane
    ring_positions_planar = ring_positions_planar[:, :-1]
    
    bond_plane_intersections @= H.T # apply alignment reflection
    assert np.allclose(bond_plane_intersections[:, -1], 0.0) # double-check that the above manipulations have aligned the array into the XY-plane
    bond_plane_intersections = bond_plane_intersections[:, :-1]

    # check for intersections inside triangulation
    tri = Delaunay(ring_positions_planar)
    pierces_ring : np.ndarray[bool] = (tri.find_simplex(bond_plane_intersections) != -1) # points not in a simplex (i.e. not intersecting the interior of the ring) are assigned to simplex "-1"
    
    return tuple(
        tuple(piercing_idx_pair.tolist()) # tolist() needed to convert from numpy int types to Python ones
            for piercing_idx_pair in spanning_pair_idxs[pierces_ring & t_is_supported] # take piercing pair to be those which intersect the ring polygon AND have interpolation parameter 0 <= t <= 1
    ) 
ring_piercing_idxs = PINPRICS # alias for convenience
    
    
# CONVENIENCE METHODS SPECIFIC TO RDKIT
def detect_ring_is_pierced(mol : Chem.Mol, ring_atom_idxs : list[int], conformer_idx : int=0) -> tuple[tuple[int, int]]:
    '''
    Detects and returns all bonded pairs of atoms which pierce a particular ring in a molecule
    
    Parameters
    ----------
    mol : Mol
        An RDKit Mol instance assumed to have at least one conformer
    ring_atom_idxs : list[int]
        List of indices of the atoms which make up the ring to be tested for piercing
    conformer_idx : int, default 0
        The index of the conformer to draw coordinates from
        By default, will just take the first conformer

    Returns
    -------
    piercing_idxs : tuple[tuple[int, int]]
        A tuple of all pairs of atom indices which were perceived to pierce the ring
    '''
    return ring_piercing_idxs(
        positions=mol.GetConformer(conformer_idx).GetPositions(), # NOTE: no need to test for presence of conformer, RDKit will already raise Exception
        bonded_atom_idxs=[(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in mol.GetBonds()], # compile array-based representation of molecule graph topology,
        ring_atom_idxs=ring_atom_idxs,
    )
    
def summarize_ring_piercing(mol : Chem.Mol, conformer_idx : int=0) -> dict[tuple[int], tuple[tuple[int, int]]]:
    '''
    Apply PINPRICS algorithm to all ring detected in a molecule and 
    provide a summary of which bonds (if any) were detected to be piercing each ring
    
    Parameters
    ----------
    mol : Mol
        An RDKit Mol instance assumed to have at least one conformer
    conformer_idx : int, default 0
        The index of the conformer to draw coordinates from
        By default, will just take the first conformer

    Returns
    -------
    piercing_summary : dict[tuple[int], tuple[tuple[int, int]]]
        A dict mapping the indices of the atoms in each ring to 
        a tuple of all pairs of atom indices which were perceived to pierce the ring
    '''
    return {
        ring_atom_idxs : detect_ring_is_pierced(mol, ring_atom_idxs, conformer_idx=conformer_idx)
            for ring_atom_idxs in mol.GetRingInfo().AtomRings()
    }