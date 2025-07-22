'''Behavioral tests for chemical reaction execution and mapping of atoms and bonds thru reactions'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import pytest

from polymerist.rdutils.reactions.reactions import AnnotatedReaction
from . import RXNS


# DEFINING REACTANTS
@pytest.fixture(scope='module')
def test_reactant_smiles() -> dict[str, str]:
    return {
        'polyester' : ('OCCO', 'O(C=O)c1ccc(cc1)C(=O)O'), # PET,
        'polyamide' :('NCCCCCCN', 'O=C(O)CCCCC(=O)O'), # Nylon-6,6
        'polyimide' : ('O(c1ccc(N)cc1)c2ccc(cc2)N', 'C1=C2C(=CC3=C1C(=O)OC3=O)C(=O)OC2=O'), # DuPont Kapton (poly (4,4'-oxydiphenylene-pyromellitimide))
        'polycarbonate_phosgene' : ('Oc1ccc(cc1)C(c2ccc(O)cc2)(C)C', 'ClC(=O)Cl'), # BPA + phosgene
        'polycarbonate_nonphosgene' : ('Oc1ccc(cc1)C(c2ccc(O)cc2)(C)C', 'O=C(Oc1ccccc1)Oc2ccccc2'), # BPA + diphenyl carbonate
        'polyurethane_isocyanate' : ('O=C=N\CCCCCC/N=C=O', 'OCCCCO'), # Bayer HDI + BDO
        'polyurethane_nonisocyanate' : ('CC(=C)C(=O)OCC1COC(=O)O1', 'NCCCCCCN'), # PCA (propylene carbonate acrylate) + hexamethylenediamine
        'polyvinyl_head_tail' : ('c1ccccc1C=C', 'c1ccccc1C=C'), # polystyrene
    }
    
