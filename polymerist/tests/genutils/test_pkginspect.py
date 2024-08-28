'''Unit tests for package inspection utilities'''

from types import ModuleType

import pytest
import math, json # use these as test cases, since they are pretty stable in stdlib

from polymerist import polymerist # this is a dummy toplevel module, and NOt the entire polymerist package
from polymerist import genutils
from polymerist.genutils import pkginspect
from polymerist.tests import data as test_data



# TABULATED EXPECTED TESTS OUTPUTS
are_modules = [
    ('--not_a_module--', False), # deliberately weird to ensure this never accidentally clashes with a legit module name
    (math, True),
    ('math', True), # test that the string -> module resolver also works as intended
    (json, True),
    ('json', True),
    (json.decoder, True),
    ('json.decoder', True),
    (polymerist, True),
    ('polymerist.polymerist', True),
    (genutils, True),
    ('polymerist.genutils', True),
]

are_packages = [
    ('--not_a_package--', False), # deliberately weird to ensure this never accidentally clashes with a legit module name
    (math, False),
    ('math', False), # test that the string -> module resolver also works as intended
    (json, True),
    ('json', True),
    (json.decoder, False),
    ('json.decoder', False),
    (polymerist, False),
    ('polymerist.polymerist', False),
    (genutils, True),
    ('polymerist.genutils', True),
]


# UNIT TESTS
@pytest.mark.parametrize('module, expected_output', are_modules)
def test_is_module(module : ModuleType, expected_output : bool) -> None:
    '''See if Python module perception behaves as expected'''
    assert pkginspect.is_module(module) == expected_output

@pytest.mark.parametrize('module, expected_output', are_packages)
def test_is_package(module : ModuleType, expected_output : bool) -> None:
    '''See if Python package perception behaves as expected'''
    assert pkginspect.is_package(module) == expected_output

