'''Unit tests for `pkginspect` package`'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from types import ModuleType

import pytest
from pathlib import Path
import math, json # use these as test cases, since they are pretty stable in stdlib

from polymerist import polymerist # this is a dummy toplevel module, and NOt the entire polymerist package
from polymerist import genutils
from polymerist.genutils import importutils
from polymerist.genutils.importutils import pkginspect
from polymerist import tests



# TABULATED EXPECTED TESTS OUTPUTS
non_module_types = [ # types that are obviously not modules OR packages, and which should fail
    bool, int, float, complex, tuple, list, dict, set, # str, Path # str and Path need to be tested separately
]

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


# MODULE AND PACKAGE PERCEPTION
@pytest.mark.parametrize('module, expected_output', are_modules)
def test_is_module(module : ModuleType, expected_output : bool) -> None:
    '''See if Python module perception behaves as expected'''
    assert pkginspect.is_module(module) == expected_output

@pytest.mark.parametrize('non_module_type', non_module_types)
def test_is_module_fail_on_invalid_types(non_module_type : type) -> None:
    '''check that module perception fails on invalid inputs'''
    with pytest.raises(AttributeError) as err_info:
        instance = non_module_type() # create a default instance
        _ = pkginspect.is_module(instance)

@pytest.mark.parametrize('module, expected_output', are_packages)
def test_is_package(module : ModuleType, expected_output : bool) -> None:
    '''See if Python package perception behaves as expected'''
    assert pkginspect.is_package(module) == expected_output

@pytest.mark.parametrize('non_module_type', non_module_types) # NOTE: these args are in fact deliberately NOT renamed to ".*package" from ".*module"
def test_is_module_fail_on_invalid_types(non_module_type : type) -> None:
    '''check that package perception fails on invalid inputs'''
    with pytest.raises(AttributeError) as err_info:
        instance = non_module_type() # create a default instance
        _ = pkginspect.is_package(instance)

# FETCHING DATA FROM PACKAGES
@pytest.mark.parametrize(
    'rel_path, module',
    [
        ('data', tests),
        ('data/sample.dat', tests),
        pytest.param('daata/simple.dat', tests, marks=pytest.mark.xfail(raises=ValueError, reason="This isn't a real file", strict=True)),
        ('pkginspect.py', importutils),
        pytest.param('fake/whatever.txt', pkginspect, marks=pytest.mark.xfail(raises=TypeError, reason="Module is not a package and therefore cannot contain resources", strict=True)),
    ]
)
def test_get_resource_path(rel_path : str, module : ModuleType) -> None:
    '''Test fetching a resource (i.e. file OR dir) from a package'''
    resource_path = pkginspect.get_resource_path_within_package(rel_path, module)
    assert isinstance(resource_path, Path)

@pytest.mark.parametrize(
    'rel_path, module',
    [
        pytest.param('data', tests, marks=pytest.mark.xfail(raises=FileNotFoundError, reason="This is a directory, NOT a file", strict=True)),
        ('data/sample.dat', tests),
        pytest.param('daata/simple.dat', tests, marks=pytest.mark.xfail(raises=ValueError, reason="This isn't a real file", strict=True)),
        ('pkginspect.py', importutils),
        pytest.param('fake/whatever.txt', pkginspect, marks=pytest.mark.xfail(raises=TypeError, reason="Module is not a package and therefore cannot contain resources", strict=True)),
    ]
)
def test_get_file_path(rel_path : str, module : ModuleType) -> None:
    '''Test fetching a file (i.e. NOT a dir) from a package'''
    resource_path = pkginspect.get_file_path_within_package(rel_path, module)
    assert isinstance(resource_path, Path)

@pytest.mark.parametrize(
    'rel_path, module',
    [
        ('data', tests),
        pytest.param('data/sample.dat', tests, marks=pytest.mark.xfail(raises=NotADirectoryError, reason='This IS a real file, but not a directory', strict=True)),
        pytest.param('daata/simple.dat', tests, marks=pytest.mark.xfail(raises=ValueError, reason="This isn't a real file", strict=True)),
        pytest.param('pkginspect.py', importutils, marks=pytest.mark.xfail(raises=NotADirectoryError, reason='This IS a real file, but not a directory', strict=True)),
        pytest.param('fake/whatever.txt', pkginspect, marks=pytest.mark.xfail(raises=TypeError, reason="Module is not a package and therefore cannot contain resources", strict=True)),
    ]
)
def test_get_dir_path(rel_path : str, module : ModuleType) -> None:
    '''Test fetching a dir (i.e. NOT a file) from a package'''
    resource_path = pkginspect.get_dir_path_within_package(rel_path, module)
    assert isinstance(resource_path, Path)