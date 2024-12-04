'''Unit tests for `attrs` package'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import pytest
from pathlib import Path

from polymerist.genutils.importutils.pkginspect import get_file_path_within_package
from polymerist.tests import data as testdata

from polymerist.polymers import building

@pytest.fixture
def fragments_path() -> Path:
    return get_file_path_within_package('peg=pla-pga.json', testdata)

# Also add separate tests module for polymers.estimation
