'''Unit tests for special SMARTS queries'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import pytest

from dataclasses import dataclass
from itertools import chain

from rdkit import Chem
from polymerist.smileslib.special import (
    _special_queries,
    SPECIAL_QUERY_ATOMS,
    SPECIAL_QUERY_MOLS,
    SPECIAL_QUERY_SMARTS,
)

# DEFINE MOLECULES TO BE QUERIED ACROSS TESTS, ALONG WITH EXPECTED QUERY RESULTS
grignard_reagent = Chem.MolFromSmiles('[Cl:1]-[Mg:2]-[C:3](-[H:4])(-[H:5])(-[H:6])', sanitize=False)
Chem.SanitizeMol(grignard_reagent) # test perfect query example here; is small but contains carbon, hydrogen, metal, halogen, and hydrogens

parabromophenyllithium = Chem.MolFromSmiles('[Br:1]-[C:2]1=[C:3](-[H:9])-[C:4](-[H:10])=[C:5](-[Li:6])-[C:7](-[H:11])=[C:8]-1-[H:12]', sanitize=False)
Chem.SanitizeMol(parabromophenyllithium) # test perfect query example here; is small but contains carbon, hydrogen, metal, halogen, and hydrogens

SPECIAL_QUERY_TEST_EXAMPLES : dict[Chem.Mol, set[int]] = {
    grignard_reagent : {
        'MH' : {2, 4, 5, 6},
        'XH' : {1, 4, 5, 6},
        'AH' : {1, 2, 3, 4, 5, 6},
        'QH' : {1, 2, 4, 5, 6},
        'M'  : {2},
        'X'  : {1},
        'A'  : {1, 2, 3},
        'Q'  : {1, 2},
    },
    parabromophenyllithium : {
        'MH' : {6, 9, 10, 11, 12},
        'XH' : {1, 9, 10, 11, 12},
        'AH' : {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12},
        'QH' : {1, 6, 9, 10, 11, 12},
        'M'  : {6},
        'X'  : {1},
        'A'  : {1, 2, 3, 4, 5, 6, 7, 8},
        'Q'  : {1, 6},
    }
}

# CLASS FOR BLUK UNIT TESTS
@pytest.mark.parametrize(
    'mol_to_query,query_name,expected_matched_map_nums',
    [
        (test_mol, query_alias, expected_matched_map_nums)
            for test_mol, expected_match_dict in SPECIAL_QUERY_TEST_EXAMPLES.items()
                for query_identifier, expected_matched_map_nums in expected_match_dict.items()
                    for query_alias in _special_queries[query_identifier]
    ]
)
class TestSpecialQueries:
    def test_special_query_atoms(self, mol_to_query : Chem.Mol, query_name : str, expected_matched_map_nums : set[int]) -> None:
        '''Test that the defined special QueryAtom definitions match the queries expected'''
        matched_map_nums = set(
            atom.GetAtomMapNum()
                for atom in mol_to_query.GetAtomsMatchingQuery(SPECIAL_QUERY_ATOMS[query_name])
        )
        assert matched_map_nums == expected_matched_map_nums
        
    def test_special_query_mols(self, mol_to_query : Chem.Mol, query_name : str, expected_matched_map_nums : set[int]) -> None:
        '''Test that the defined special Chem.Mol definitions match the queries expected'''
        matched_map_nums = set(
            mol_to_query.GetAtomWithIdx(atom_idx).GetAtomMapNum()
                for atom_idx in chain(*mol_to_query.GetSubstructMatches(SPECIAL_QUERY_MOLS[query_name]))
        )
        assert matched_map_nums == expected_matched_map_nums
        
    def test_special_query_smarts(self, mol_to_query : Chem.Mol, query_name : str, expected_matched_map_nums : set[int]) -> None:
        '''Test that the defined special SMARTS definitions match the queries expected'''
        query_mol = Chem.MolFromSmarts(SPECIAL_QUERY_SMARTS[query_name])
        matched_map_nums = set(
            mol_to_query.GetAtomWithIdx(atom_idx).GetAtomMapNum()
                for atom_idx in chain(*mol_to_query.GetSubstructMatches(query_mol))
        )
        assert matched_map_nums == expected_matched_map_nums