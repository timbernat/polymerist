'''Unit tests for `dependencies` package'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import pytest

from typing import Any, Callable
from polymerist.genutils.importutils import dependencies


# Testing module finding
@pytest.mark.parametrize(
    'module_names, expected_found', [
        (['polymerist'], True),     # we'd better hope the parent module is present if we're running tests on it :P
        (['sys'], True),            # test stdlib packages which ought to be present if Python is
        (['os', 'sys'], True),      # test that unpacking also works
        (['fake--module'], False),  # test an obviously fake module name (don't want to try an actual module in case it becomes an dependency someday)
        ([42], False),              # test something that isn't even a module to check error handling
    ]
)
def test_modules_installed(module_names : list[str], expected_found : bool) -> None:
    '''Check that module install checker correctly identifies present and absent modules'''
    assert dependencies.modules_installed(*module_names) == expected_found
    
# Testing requires_modules decorator
@dependencies.requires_modules('os')
def should_pass() -> str:
    '''Dummy function to test requires_modules decorator for dependencies that are present'''
    return 'I will run!'

@dependencies.requires_modules('fake--module')
def should_fail() -> str:
    '''Dummy function to test requires_modules decorator for dependencies that are present'''
    return 'I will xfail :('

@pytest.mark.parametrize(
    'func',
    [
        should_pass,
        pytest.param(should_fail, marks=pytest.mark.xfail(raises=ImportError, reason='The required module shouldn\'t be found in the environment', strict=True)),
    ]
)
def test_requires_modules(func : Callable[..., Any]) -> None:
    '''Test that the requires_modules decortor correctly wraps functions'''
    _ = func() # no assertion needed, xfail cases should raise Exception while working cases will ternimate without Exception