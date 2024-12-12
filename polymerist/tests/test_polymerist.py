"""
Unit and regression test for the polymerist package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import polymerist


def test_polymerist_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "polymerist" in sys.modules
