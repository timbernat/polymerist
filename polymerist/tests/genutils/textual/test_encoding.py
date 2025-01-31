'''Unit tests for `encoding` package'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import pytest

from polymerist.genutils.textual.encoding import representable_as_int


@pytest.mark.parametrize(
    'string',
    [
        # test signed and unsigned 0
        '0',
        '+0',
        '-0',
        # test single-digit integers
        '1',
        '5',
        '-9',
        '+6',
        # test multi-digit integers
        '12',
        '1232164',
        '+192',
        '-1232',
        # test integers with underscores (Python exclusive)
        '125_000',
        '1_2_3_4',
        '+1_00_00_000', # one "crore"
        '-3_14_1592_6',
    ],
)
def test_int_repr_valid(string : str) -> None:
    '''Test that string which represent valid integers are correctly identified'''
    assert representable_as_int(string) == True

@pytest.mark.parametrize(
    'string',
    [
        # test rejection of leading zeros for integers not identically 0
        '01',
        '00004241',
        '010101',
        '+0071',
        '+01_618',
        '-0002',
        # test rejection on string with non-numeric characters
        '', # in particular, the empty string should not be recognized
        'fbhf',
        'a123',
        'fubar',
        'fubar69',
        '+a5sd',
        '-123x232',
        '2+2',
        # test rejection of ints with leading or tailing underscores
        '-10_',
        '_99',
        '+_0139_',
    ],
)
def test_int_repr_invalid(string : str) -> None:
    '''Test that string which represent invalid integers are correctly rejected'''
    assert representable_as_int(string) == False