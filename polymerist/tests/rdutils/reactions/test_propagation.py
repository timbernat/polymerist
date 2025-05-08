'''Behavioral tests to ensure reaction fragment enumeration behaves as expected'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import pytest

from polymerist.rdutils.reactions.reactions import AnnotatedReaction
from polymerist.rdutils.reactions.reactors import PolymerizationReactor
from polymerist.rdutils.reactions.fragment import ReseparateRGroups, CutMinimumCostBondsStrategy
from . import RXNS


# DEFINING REACTANTS
@pytest.fixture(scope='module')
def test_reactant_smiles() -> dict[str, str]:
    # TOSELF: Ought to include at least one example of the many logical edge cases I'm trying to
    # handle (e.g. autopolymerization, large monomers, more than 2 monomers, etc.)
    multifunct_examples : list[str] = [
        'NCCCCCC(=O)O.NCCCCCCN.O=C(O)CCCCC(=O)O',
        'NCCCCCC(=O)O.NCCCCCCN.O=C(O)c1ccc(C(=O)O)cc1',
        'CC(C)(c1ccc(O)cc1)c1ccc(O)cc1.O=C(Cl)Cl.Oc1ccc(Oc2ccc(O)cc2)cc1',
        'O=C(O)CCC(=O)O.O=C(O)c1ccc(C(=O)O)o1.OCCCCO',
    ]
    ...