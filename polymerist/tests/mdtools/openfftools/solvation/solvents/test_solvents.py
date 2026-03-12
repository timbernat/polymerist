'''Testing that pre-compiled solvents are defined correctly'''

import pytest

from numpy.testing import assert_equal
from openff.toolkit import Molecule
from polymerist.mdtools.openfftools.solvation.solvents import water_TIP3P

@pytest.mark.parametrize(
        'solvent,expected_charge',
        [
            (water_TIP3P, 0.0),
            # DEV: left as parametrize in case of addition of other solvents in the future
        ]
)
def test_solvent_net_charge(solvent : Molecule, expected_charge : float) -> None:
    '''Test that solvents actually have the expected net charge'''
    assert_equal(solvent.partial_charges.sum().m_as('elementary_charge'), expected_charge)