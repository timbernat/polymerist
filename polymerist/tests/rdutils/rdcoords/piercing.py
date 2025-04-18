'''Unit tests for `sanitization` package'''

import pytest
from pathlib import Path

from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdmolops import SanitizeMol
from rdkit.Chem.rdmolfiles import SDMolSupplier

from polymerist.genutils.importutils.pkginspect import get_dir_path_within_package
from polymerist.tests import data as testdata

from polymerist.rdutils.rdcoords.piercing import detect_ring_is_pierced, summarize_ring_piercing


@pytest.fixture
def test_structures_dir_path() -> Path:
    return get_dir_path_within_package('ring_piercing_examples', testdata)

@pytest.fixture
def unpierced_structure(test_structures_dir_path : Path) -> Mol:
    '''A control structure which has NO ring piercings'''
    with SDMolSupplier(str(test_structures_dir_path / 'unpierced.sdf'), sanitize=False, removeHs=False) as suppl:
        unpierced_mol = suppl[0]
        SanitizeMol(unpierced_mol)
        
    return unpierced_mol

@pytest.fixture
def pierced_structure(test_structures_dir_path : Path) -> Mol:
    '''A test structure which has many interlocked ring piercings'''
    with SDMolSupplier(str(test_structures_dir_path / 'pierced.sdf'), sanitize=False, removeHs=False) as suppl:
        pierced_mol = suppl[0]
        SanitizeMol(pierced_mol)
        
    return pierced_mol


def test_control(unpierced_structure : Mol) -> None:
    '''Test that structures with NO ring piercings are not misidentified as having piercing events'''
    piercing_map = summarize_ring_piercing(unpierced_structure)
    assert not any(piercing_map.values()) # check that indices of all piercing bonds are empty
    
@pytest.mark.parametrize(
    'ring_atom_idxs,expected_piercing_idxs',
    [
        [ (94, 95, 97, 99, 101, 103), set() ], # use frozensets to allow "set of sets" (desirable, since only the CONTENT, not the order of piercing indices matters
        [ (71, 72, 74, 76, 78, 80)  , set() ],
        [ (10, 9, 18, 16, 14, 12)   , set({frozenset({36, 38})}) ],
        [ (62, 53, 54, 56, 58, 60)  , set({frozenset({31, 40})}) ],
        [ (38, 40, 31, 32, 34, 36)  , set({frozenset({18, 9}), frozenset({58, 59}), frozenset({58, 56})}) ],
    ]
)
def test_piercing_detection(pierced_structure : Mol, ring_atom_idxs : tuple[int], expected_piercing_idxs : set[frozenset[int, int]]) -> None:
    '''Test that the ring piercing algorithm implemented correctly detects 
    all instances and locations of piercing in the example pierced structure'''
    pierced_idxs_found = set(
        frozenset(idx_pair) # convert to sets to enable order-independent equality check
            for idx_pair in detect_ring_is_pierced(pierced_structure, ring_atom_idxs, conformer_idx=0)
    )
    assert expected_piercing_idxs == pierced_idxs_found
