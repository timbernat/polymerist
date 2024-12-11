'''Unit tests for `substrings` package'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import pytest

from polymerist.genutils.textual.substrings import unique_string, shortest_repeating_substring, repeat_string_to_length


@pytest.mark.parametrize('string', ['Lorem', 'ipsum', 'dolor', 'sit', 'amet', 'consectetur', 'adipiscing', 'elit'])
def test_unique_str_unordered(string : str) -> None:
    '''Test that unique characters are coorectly identified WITHOUT respect to order'''
    assert set(unique_string(string, preserve_order=False)) == set(string)

@pytest.mark.parametrize('string, expected_output',
    [
        ('aaaaa', 'a'),
        ('BABAB', 'BA'),
        ('balaclava', 'balcv'),
        ('catamaran', 'catmrn'),
        ('unique', 'uniqe'), # self-reference makes everything better :P
        ('singular', 'singular'), # test string with aready-unique characters are unaffected
    ]
)
def test_unique_str_ordered(string : str, expected_output : str) -> None:
    '''Test that unique characters are coorectly identified WITH respect to order'''
    assert unique_string(string, preserve_order=True) == expected_output
    
    
@pytest.mark.parametrize('string, expected_output',
    [
        ('aaaaa', 'a'),
        ('booboo', 'boo'),
        ('piripiri', 'piri'),
        ('ababab', 'ab'),
        ('bcbabcbabcba', 'bcba'),
        # sequences which do not repeat a whole-number of times
        ('no repeats', 'no repeats'),
        ('ababa', 'ababa'), 
        ('bonobo', 'bonobo'),
    ]
)
def test_shortest_repeating_substring(string : str, expected_output : str) -> None:
    '''Test that minimal repeating substrings are correctly identified'''
    assert shortest_repeating_substring(string) == expected_output
    
    
@pytest.mark.parametrize('string, target_length, joiner, expected_output',
    [
        ('BACA', 10, '', 'BACABACABA'), # expected "standard" use case
        ('BACA', 1, '', 'B'),           # test case where target length is shorter than the whole string
        ('BACA', 0, '', ''),            # test that no repeats yields the empty string
        ('BACA', 4, '', 'BACA'),             # test precisely one repeat without joins
        ('BACA', 10, '|', 'BACA|BACA|BA'),   # test joiners
        ('BACA', 4, '|', 'BACA'),            # test no joiners are added when exactly one string repeat occurs
        ('BACA', 12, '|', 'BACA|BACA|BACA'), # test no extraneous joiners are included for purely-whole number of repeats
        ('CAT', 5, '', 'CATCA'), # test with triads (and different base string)
        pytest.param(''   ,   7, '', None, marks=pytest.mark.xfail(raises=ValueError, reason='Empty string can\'t be repeated into nonempty string', strict=True)),
        pytest.param('CAT', 4.2, '', None, marks=pytest.mark.xfail(raises=TypeError , reason='Non-integer string length doesn\'t make sense', strict=True)),
        pytest.param('CAT',  -1, '', None, marks=pytest.mark.xfail(raises=IndexError, reason='Can\'t have string with fewer than 0 characters', strict=True)),
    ]
)
def test_repeat_string_to_length(string : str, target_length : int, joiner : str, expected_output : str) -> None:
    '''Test that string repetition to a given length returns the expected string WITH joingin characters present'''
    assert repeat_string_to_length(string, target_length=target_length, joiner=joiner) == expected_output