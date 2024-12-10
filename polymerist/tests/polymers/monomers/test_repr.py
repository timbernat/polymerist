'''Tests that collections of monomer fragments are treated as expected'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Any

import pytest

from ..import PE_FRAGMENTS, MPD_TMC_FRAGMENTS, PEG_PLGA_FRAGMENTS
from polymerist.polymers.monomers.repr import MonomerGroup


# Example fragments groups
@pytest.fixture(scope='function') # want to re-initialize for each test function to avoid cross-contamination
def monogrp_degenerate() ->  MonomerGroup:
    return MonomerGroup(monomers={})

@pytest.fixture(scope='function')
def monogrp_polyethylene() ->  MonomerGroup:
    return MonomerGroup(monomers=PE_FRAGMENTS)

@pytest.fixture(scope='function')
def monogrp_mpd_tmc() ->  MonomerGroup:
    return MonomerGroup(monomers=MPD_TMC_FRAGMENTS)

@pytest.fixture(scope='function')
def monogrp_peg_plga() ->  MonomerGroup:
    return MonomerGroup(monomers=PEG_PLGA_FRAGMENTS)


# Testing all routes to initialization
@pytest.mark.parametrize(
    'monomers',
    [
        { # nominal test case
            'PGA-1A': ['[#8D2+0:1](-[#6D4+0:2](-[#6D3+0:3](=[#8D1+0:4])-[*:5])(-[#1D1+0:7])-[#1D1+0:8])-[#1D1+0:6]'],
            'PGA-1B': ['[*:1]-[#8D2+0:2]-[#6D4+0:3](-[#6D3+0:4](=[#8D1+0:5])-[#8D2+0:6]-[#1D1+0:9])(-[#1D1+0:7])-[#1D1+0:8]'],
            'PGA-2': ['[*:1]-[#8D2+0:2]-[#6D4+0:3](-[#6D3+0:4](=[#8D1+0:5])-[*:6])(-[#1D1+0:7])-[#1D1+0:8]'],
        },
        { # test that list closure autofill works
            'PGA-1A': '[#8D2+0:1](-[#6D4+0:2](-[#6D3+0:3](=[#8D1+0:4])-[*:5])(-[#1D1+0:7])-[#1D1+0:8])-[#1D1+0:6]',
            'PGA-1B': '[*:1]-[#8D2+0:2]-[#6D4+0:3](-[#6D3+0:4](=[#8D1+0:5])-[#8D2+0:6]-[#1D1+0:9])(-[#1D1+0:7])-[#1D1+0:8]',
            'PGA-2': '[*:1]-[#8D2+0:2]-[#6D4+0:3](-[#6D3+0:4](=[#8D1+0:5])-[*:6])(-[#1D1+0:7])-[#1D1+0:8]',
        },
        # XFAILS: test that the initializer rejects...
        pytest.param(
            { # ...1) non-string like objects
                'foo' : 42.0,
                'bar' : True,
            },
            marks=pytest.mark.xfail(
                raises=TypeError,
                reason='Monomer fragment inputs are not stringlike',
                strict=True,
            ),
        ),
        pytest.param(
            { #  1a) more subtly, list OF CONTAINERS of valid SMARTS are still invalid
                'PGA-1A': [( 
                    '[#8D2+0:1](-[#6D4+0:2](-[#6D3+0:3](=[#8D1+0:4])-[*:5])(-[#1D1+0:7])-[#1D1+0:8])-[#1D1+0:6]',
                    '[#8D2+0:1](-[#6D4+0:2](-[#6D3+0:3](=[#8D1+0:4])-[*:5])(-[#1D1+0:7])-[#1D1+0:8])-[#1D1+0:6]'
                )],
                'PGA-2': ['[*:1]-[#8D2+0:2]-[#6D4+0:3](-[#6D3+0:4](=[#8D1+0:5])-[*:6])(-[#1D1+0:7])-[#1D1+0:8]'],
            },
            marks=pytest.mark.xfail(
                raises=TypeError,
                reason='Monomer fragment inputs are not stringlike',
                strict=True,
            ),
        ),
        pytest.param(
            { # ...2) empty lists
                'PGA-1A': [],
                'PGA-2' : [],
            },
            marks=pytest.mark.xfail(
                raises=IndexError,
                reason='At least one monomer fragment input is empty',
                strict=True,
            ),
        ),
        pytest.param(
            { # ...3) non-empty strings which are nevertheless invalid SMARTS 
             #- NOTE: empty strings, perhaps surprisingly, actually ARE valid as SMARTS and therefore aren't xfail tested here
                'fake-1': ['this is a bogus SMARTS'],
                'invalid-2': ['so_is_this'],
            },
            marks=pytest.mark.xfail(
                raises=ValueError,
                reason='At least one monomer fragment input is not valid a SMARTS string',
                strict=True,
            ),
        ),
        pytest.param(
            { # ...3a) this one is very subtle, but SMARTS with slight errors which invalidate themas SMARTS should also fail
                'PGA-1A': ['[OH]CD(=O)*'], # fat-finger mistake, "D" should be "C"
                'PGA-2': ['*OCC(+O)*'], # forgot to hit shift when typing double bond
            },
            marks=pytest.mark.xfail(
                raises=ValueError,
                reason='At least one monomer fragment input is not valid a SMARTS string',
                strict=True,
            ),
        ),
    ]
)
def test_monogrp_init(monomers : dict[str, Any]) -> None:
    '''Check that the MonomerGroup initializer handles valid inputs as expected and reject invalid inputs'''
    _ = MonomerGroup(monomers=monomers) # no assert needed, just checking when initialization completes
    
# Testing properties of contained monomers
@pytest.mark.parametrize(
    'monogrp, expected_is_linear',
    [
        ('monogrp_peg_plga', True),
        ('monogrp_mpd_tmc', False),
    ],
)
def test_monogrp_linearity(monogrp : MonomerGroup, expected_is_linear : bool, request : pytest.FixtureRequest) -> None:
    '''Test whether branched and unbranched chain fragment detection behaves as expected'''
    monogrp = request.getfixturevalue(monogrp) # unpack fixtures into their respective values
    assert monogrp.is_linear == expected_is_linear
    
@pytest.mark.parametrize(
    'monogrp, expected_counts',
    [
        ('monogrp_peg_plga', (3, 6)),
        ('monogrp_mpd_tmc', (3, 2)),
    ],
)
def test_monogrp_mid_and_term_counts(monogrp : MonomerGroup, expected_counts : tuple[int, int], request : pytest.FixtureRequest) -> None:
    '''Test whether middle and terminal monomers are counted correctly'''
    monogrp = request.getfixturevalue(monogrp) # unpack fixtures into their respective values
    assert monogrp.num_mid_and_term == expected_counts
    
# Testing end group determination
@pytest.mark.parametrize(
    'monogrp, term_orient, expected_end_groups',
    [
        ('monogrp_peg_plga', {}, {'head' : 'PEG-1A', 'tail' : 'PEG-1B'}), # test autogeneration from first 2 when
        ('monogrp_mpd_tmc', (3, 2)),
    ],
)
def test_monogrp_end_groups(monogrp : MonomerGroup, term_orient : dict[str, str], expected_end_groups : dict[str, str], request : pytest.FixtureRequest) -> None:
    '''Test whether procedural end group determination'''
    monogrp = request.getfixturevalue(monogrp) # unpack fixtures into their respective values
    monogrp.term_orient = term_orient
    
    end_group_catalogue = monogrp.linear_end_groups()
    end_group_names = {
        head_or_tail : resname  # drop RDKit Mol for check (Mol object is harder to validate, use bound name as proxy)
            for head_or_tail, (resname, _) in end_group_catalogue.items()
    }
    
    assert end_group_names == expected_end_groups
    