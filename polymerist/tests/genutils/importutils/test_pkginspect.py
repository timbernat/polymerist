'''Unit tests for `pkginspect` package`'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import pytest
from _pytest.mark.structures import ParameterSet

from types import ModuleType
from dataclasses import dataclass

from pathlib import Path
import math, json # use these as test cases, since they are pretty stable in stdlib

from polymerist import polymerist # this is a dummy toplevel module, and NOt the entire polymerist package
from polymerist import genutils
from polymerist.genutils import importutils
from polymerist.genutils.importutils import pkginspect
from polymerist import tests


# TABULATED EXPECTED TESTS OUTPUTS
def non_module_types() ->list[type]:
    '''Types that are obviously not modules OR packages, and which should fail'''
    return [
        bool, int, float, complex, tuple, list, dict, set, 
        # str, Path # str and Path need to be tested separately
    ]

def obviously_fake_resource_param() ->ParameterSet:
    return pytest.param(
        'fake/whatever.txt', pkginspect,
        marks=pytest.mark.xfail(
            raises=(TypeError, ValueError), # DEVNOTE: annoyingly, Exception raised is TypeError in Python 3.11 but ValueError in 3.12
            reason="Module is not a package and therefore cannot contain resources",
            strict=True,
        )
    )

@dataclass(frozen=True)
class ModuleExample:
    '''For encapsulating package and module check tests'''
    resource : str | ModuleType
    is_module : bool
    is_package : bool

def module_examples() -> tuple[ModuleExample, ...]:
    '''
    Module-like objects labelled with whether they
    are modules and whether they are packages
    '''
    return (
        # deliberately weird to ensure this never accidentally clashes with a legit module name
        ModuleExample('--not_a_module--', False, False), 
        ModuleExample(math, True, False),
        # test that the string -> module resolver also works as intended
        ModuleExample('math', True, False), 
        ModuleExample(json, True, True),
        ModuleExample('json', True, True),
        ModuleExample(json.decoder, True, False),
        ModuleExample('json.decoder', True, False),
        ModuleExample(polymerist, True, False),
        ModuleExample('polymerist.polymerist', True, False),
        ModuleExample(genutils, True, True),
        ModuleExample('polymerist.genutils', True, True),
    )


# MODULE AND PACKAGE PERCEPTION
@pytest.mark.parametrize('non_module_type', non_module_types())
def test_invalid_module_type_rejected(non_module_type : type) -> None:
    '''Check that module loading on obviously non-module types raises explicit Exception'''
    with pytest.raises(ModuleNotFoundError) as err_info:
        instance = non_module_type() # create a default instance
        _ = pkginspect._load_module(instance)

@pytest.mark.parametrize('module_example', module_examples())
def test_is_module(module_example : ModuleExample) -> None:
    '''See if Python module perception behaves as expected'''
    assert pkginspect.is_module(module_example.resource) == module_example.is_module

@pytest.mark.parametrize('non_module_type', non_module_types())
def test_is_module_fail_on_invalid_types(non_module_type : type) -> None:
    '''Check that module perception fails on invalid input types'''
    instance = non_module_type() # create a default instance
    assert not pkginspect.is_module(instance)

@pytest.mark.parametrize('module_example', module_examples())
def test_is_package(module_example : ModuleExample) -> None:
    '''See if Python package perception behaves as expected'''
    assert pkginspect.is_package(module_example.resource) == module_example.is_package

@pytest.mark.parametrize('non_package_type', non_module_types())
def test_is_package_fail_on_invalid_types(non_package_type : type) -> None:
    '''Check that package perception fails on invalid input types'''
    instance = non_package_type() # create a default instance
    assert not pkginspect.is_package(instance)

# FETCHING DATA FROM PACKAGES
@pytest.mark.parametrize(
    'rel_path, module',
    [
        ('data', tests),
        ('data/sample.dat', tests),
        pytest.param(
            'daata/simple.dat', tests,
            marks=pytest.mark.xfail(
                raises=ValueError,
                reason="This isn't a real file",
                strict=True
            ),
        ),
        ('pkginspect.py', importutils),
        obviously_fake_resource_param(),
    ]
)
def test_get_resource_path(rel_path : str, module : ModuleType) -> None:
    '''Test fetching a resource (i.e. file OR dir) from a package'''
    resource_path = pkginspect.get_resource_path_within_package(rel_path, module)
    assert isinstance(resource_path, Path)

@pytest.mark.parametrize(
    'rel_path, module',
    [
        pytest.param(
            'data', tests,
            marks=pytest.mark.xfail(
                raises=FileNotFoundError,
                reason="This is a directory, NOT a file",
                strict=True,
            )
        ),
        ('data/sample.dat', tests),
        pytest.param(
            'daata/simple.dat', tests,
            marks=pytest.mark.xfail(
                raises=ValueError,
                reason="This isn't a real file",
                strict=True,
            )
        ),
        ('pkginspect.py', importutils),
        obviously_fake_resource_param(),
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
        pytest.param(
            'data/sample.dat', tests,
            marks=pytest.mark.xfail(
                raises=NotADirectoryError,
                reason='This IS a real file, but not a directory',
                strict=True,
            )
        ),
        pytest.param(
            'daata/simple.dat', tests,
            marks=pytest.mark.xfail(
                raises=ValueError,
                reason="This isn't a real file",
                strict=True,
            )
        ),
        pytest.param(
            'pkginspect.py', importutils, 
            marks=pytest.mark.xfail(
                raises=NotADirectoryError,
                reason='This IS a real file, but not a directory',
                strict=True,
            )
        ),
        obviously_fake_resource_param(),
    ]
)
def test_get_dir_path(rel_path : str, module : ModuleType) -> None:
    '''Test fetching a dir (i.e. NOT a file) from a package'''
    resource_path = pkginspect.get_dir_path_within_package(rel_path, module)
    assert isinstance(resource_path, Path)