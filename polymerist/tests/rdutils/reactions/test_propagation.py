'''Behavioral tests to ensure reaction fragment enumeration behaves as expected'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import pytest

from polymerist.rdutils.reactions.reactions import AnnotatedReaction
from polymerist.rdutils.reactions.reactors import PolymerizationReactor
from polymerist.rdutils.reactions.fragment import ReseparateRGroupsUnique, CutMinimumCostBondsStrategy
from . import RXNS


# DEFINING REACTANTS
@pytest.fixture(scope='module')
def test_reactant_smiles() -> dict[str, str]:
    # TOSELF: Ought to include at least one example of the many logical edge cases I'm trying to
    # handle (e.g. autopolymerization, large monomers, more than 2 monomers, etc.)
    ...