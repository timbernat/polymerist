'''Units tests from SMILES/SMARTS parsing and sanitization'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import pytest

from polymerist.smileslib.cleanup import (
    Smiles,
    Smarts,
    is_valid_SMILES,
    is_valid_SMARTS,
    expanded_SMILES,
)


@pytest.mark.parametrize(
    'smiles,is_valid_expected',
    [ # NOTE: not testing chemical validity; just SMILES-compliance
        # acetone examples
        ('CC(=O)C', True),
        ('C[C:1](=O)C', True),              # partially mapped
        ('[C:0][C:1](=[O:2])[C:3]', True),  # acetone - fully mapped
        ('[#6]C(=[#8])C', True),            # test atomic numbers (partial)
        ('[#6][#6](=[#8])[#6]', True),      # test atomic numbers (complete)
        ('[13C]C(=O)C', True),              # test isotope
        ('CC(=OC', False),                  # missing closing parenthesis
        ('CC:1](=O)C', False),              # unclosed map
        # benzoic acid examples
        ('c1ccccc1C(=O)O', True),           # benzoic acid test aromatic
        ('C1=CC=CC=C1C(=O)O', True),        # test aromatic, kekulized
        ('c1ccCcc1C(=O)O', False),          # unkekulizable
        ('c1ccccc1C(=0)O', False),          # surprise zero >:)
        # non-examples
        ('obviously-bogus', False),         # self-explanatory
        ('C6H12O6', False),
    ]
)
def test_validation_SMILES(smiles : Smiles, is_valid_expected : bool) -> None:
    '''Test that valid and invalid SMILES strings are correctly differentiated'''
    assert is_valid_SMILES(smiles) == is_valid_expected

@pytest.mark.parametrize(
    'smarts,is_valid_expected',
    [
        # nitrile examples
        ('*CN', True),
        ('*-C#N', True),
        ('[*]-[#6D2]#[#7D1]', True),
        ('[*]-[#6D2]#[#ND1]', False),   # N is not a valid atomic unmber
        ('[*]-[#6D2]#[#7D1', False),    # missing close brackets
        ('*-[nitrile]', False),
        ('obviously-bogus', False),     # self-explanatory
        # alcohol examples
        ('*-O', True),
        ('*-O-[H]', True), # make bonds/atoms explicit
        ('*-[#6X2]-[#1]', True),
        ('*-[#6X2][#1]', True),
        ('*-[#6X2-[#1]', False), # mssing bracket
        ('*-[$([#6X2]-[#1])]', True), # merge OH into single lookaround
        ('*-[([#6X2]-[#1])]', False), # lookaround without "$" delimiter is invalid
        ('[#1:0]-[#8:1]-[!$([#6]=[#8]):2]',True), # non-naive query w/ negative lookaround (i.e. check OH is not part of COOH)
        ('I\'M NOT A SMARTS!', False),
    ]
)
def test_validation_SMARTS(smarts : Smarts, is_valid_expected : bool) -> None:
    '''Test that valid and invalid SMARTS strings are correctly differentiated'''
    assert is_valid_SMARTS(smarts) == is_valid_expected
    
@pytest.mark.parametrize(
    'input_smiles,assign_map_nums,kekulize,canonicalize,smiles_expected',
    [
        # acetone
        ('CC(=O)C',False,True,True,'[H]-[C](-[H])(-[H])-[C](=[O])-[C](-[H])(-[H])-[H]'),
        ('CC(=O)C',False,True,False,'[C](-[C](=[O])-[C](-[H])(-[H])-[H])(-[H])(-[H])-[H]'), # preserve atom order (inserted Hs at end)
        ('CC(=O)C',False,False,True,'[H]-[C](-[H])(-[H])-[C](=[O])-[C](-[H])(-[H])-[H]'),   # kekulization doesn't matter here
        ('CC(=O)C',True,True,True,'[C:1](-[C:2](=[O:3])-[C:4](-[H:8])(-[H:9])-[H:10])(-[H:5])(-[H:6])-[H:7]'),  # add map numbers - mapped in order of insertion
        ('CC(=O)C',True,True,False,'[C:1](-[C:2](=[O:3])-[C:4](-[H:8])(-[H:9])-[H:10])(-[H:5])(-[H:6])-[H:7]'), # canonicalization doesn't seem to matter for mapped atoms
        # thiophene - to test kekulization
        ('s1cccc1',False,False,False,'[s]1:[c](-[H]):[c](-[H]):[c](-[H]):[c]:1-[H]'), # preserve atom order (inserted Hs at end)
        ('s1cccc1',False,True,False,'[S]1-[C](-[H])=[C](-[H])-[C](-[H])=[C]-1-[H]'),  # preserve atom order and kekulize
        ('s1cccc1',False,True,True,'[H]-[C]1=[C](-[H])-[C](-[H])=[C](-[H])-[S]-1'),   # canonicalize  + kekulize - harder to predict order
        ('s1cccc1',True,False,False,'[s:1]1:[c:2](-[H:6]):[c:3](-[H:7]):[c:4](-[H:8]):[c:5]:1-[H:9]'), # map aromatic
        ('s1cccc1',True,True,False,'[S:1]1-[C:2](-[H:6])=[C:3](-[H:7])-[C:4](-[H:8])=[C:5]-1-[H:9]'),  # map kekulized
    ]
)
def test_SMILES_expansion(
        input_smiles : Smiles,
        assign_map_nums : bool,
        kekulize : bool,
        canonicalize : bool,
        smiles_expected : Smiles,
    ) -> None:
    '''Test that "expanded" SMILES string convey the expected chemical information'''
    assert expanded_SMILES(
        input_smiles,
        assign_map_nums=assign_map_nums,
        kekulize=kekulize,
        canonicalize=canonicalize,
    ) == smiles_expected