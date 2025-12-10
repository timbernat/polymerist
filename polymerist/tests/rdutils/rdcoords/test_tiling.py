'''Unit tests for `tiling` package'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import pytest
from pathlib import Path

import numpy as np

from rdkit.Chem import Mol, SDMolSupplier, AssignStereochemistryFrom3D
from rdkit.Chem.rdMolTransforms import ComputeCanonicalTransform, TransformConformer

from polymerist.genutils.importutils.pkginspect import get_dir_path_within_package
from polymerist.tests import data as testdata


@pytest.fixture
def test_structures_dir_path() -> Path:
    return get_dir_path_within_package('stereo_inversion', testdata)

@pytest.fixture
def chiral_mol(test_structures_dir_path : Path) -> Mol: 
    # DEV: would've liked to have a SMILES-based MRE, but is difficult to know in advance when inversion will happen;
    # using example "from the wild" which has been observed to reliably exhibit this bad behavior
    with SDMolSupplier(str(test_structures_dir_path / 'trimer.sdf'), sanitize=True) as suppl:
        return suppl[0]

def test_canon_transform_stereo_inversion(chiral_mol : Mol) -> None:
    '''
    Test that RDKit's principal axis alignment doesn't inadvertently
    invert coordinate handedness and therefore stereochemistry, 
    or that at least this is caught and corrected when it happens

    Made necessary by an RDKit bug, (only patched in RDKit 2025.09.x) reported in
    https://github.com/rdkit/rdkit/issues/8992 and even earlier in https://github.com/rdkit/rdkit/issues/8720 
    ''' # TODO: invoke tiler to check that internal patches cover this issue (current test only test for underlying issue, which my code can't fix)
    chiral_atom_idx : int = 1 # atom known to be stereocenter in this structure
    conformer = chiral_mol.GetConformer(0)
    centering = ComputeCanonicalTransform(conformer)
    orient : float = np.sign(np.linalg.det(centering))

    chiral_tag_init = chiral_mol.GetAtomWithIdx(chiral_atom_idx).GetChiralTag()
    TransformConformer(conformer, centering)
    AssignStereochemistryFrom3D(chiral_mol) # ensure that stereo is updated post-transform
    chiral_tag_final = chiral_mol.GetAtomWithIdx(chiral_atom_idx).GetChiralTag()

    # DEV: not checking eps neighborhood of 1.0 to avoid precision headaches; positivity sufficies for this check 
    assert (orient > 0.0) and (chiral_tag_init == chiral_tag_final)