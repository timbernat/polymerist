'''Unit tests for filetree operations'''

import pytest
from pathlib import Path

from polymerist.genutils.pkginspect import get_dir_path_within_package
from polymerist.tests import data as testdata


@pytest.fixture
def testdir() -> Path:
    return get_dir_path_within_package('dummy_dir', testdata)
    