'''Behavioral tests for chemical reaction execution and mapping of atoms and bonds thru reactions'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Iterable
Smiles = str
Smarts = str

import pytest
from . import RXNS

from rdkit import Chem
from polymerist.rdutils.reactions.reactions import AnnotatedReaction
from polymerist.rdutils.sanitization import explicit_mol_from_SMILES, SANITIZE_ALL, AROMATICITY_MDL


def explicit_mols_from_SMILES(allsmiles : Iterable[Smiles]) -> tuple[Chem.Mol, ...]:
    '''Convert a in iterable of SMILES strings to explicit, sanitized Mols'''
    return tuple([
        explicit_mol_from_SMILES(
            smiles,
            assign_map_nums=False,
            sanitize_ops=SANITIZE_ALL,
            aromaticity_model=AROMATICITY_MDL,
        )
            for smiles in allsmiles
    ])

@pytest.mark.parametrize(
    'rxn_smarts,reactant_smiles,expected_product_patterns',
    [
        ## polyester condensation - PET
        (
            '[#1]-[#8:1]-[!$([#6]=[#8]):2].[#8](-[#6:3](=[#8:4])-[*:5])-[#1]>>[#8:1](-[!$([#6]=[#8]):2])-[#6:3](=[#8:4])-[*:5]',
            ('OCCO', 'O(C=O)c1ccc(cc1)C(=O)O'),
            ('O=COc1ccc(C(=O)OCCO)cc1',),
        ),
        ## polyamide condensation - HMDA and adipic acid into Nylon-6,6
        (
            "[#7:1](-[*:2])(-[#1])-[#1:3].[#8](-[#6:4](=[#8:5])-[*:6])-[#1]>>[#7:1](-[*:2])(-[#1:3])-[#6:4](=[#8:5])-[*:6]",
            ('NCCCCCCN', 'O=C(O)CCCCC(=O)O'),
            ('NCCCCCCNC(=O)CCCCC(=O)O',),
        ),
        ## polyimide - DuPont Kapton (poly (4,4'-oxydiphenylene-pyromellitimide))
        (
            "[#7:1](-[*:2])(-[#1])-[#1].[*:3]-[#6:4](=[#8:5])-[#8]-[#6:6](=[#8:7])-[*:8]>>[#7:1](-[*:2])(-[#6:4](-[*:3])=[#8:5])-[#6:6](=[#8:7])-[*:8]",
            ('c1cc(N)ccc1Oc1ccc(N)cc1', 'c1c2C(=O)OC(=O)c2cc3C(=O)OC(=O)c31'),
            ('c1c2C(=O)OC(=O)c2cc3C(=O)N(c4ccc(Oc5ccc(N)cc5)cc4)C(=O)c31',),
        ),
        ## polycarbonate (phosgene route) - BPA + phosgene
        (
            "[#1]-[#8:1]-[!$([#6]=[#8]):2].[#17]-[#6:3](=[#8:4])-[*:5]>>[#8:1](-[!$([#6]=[#8]):2])-[#6:3](=[#8:4])-[*:5]",
            ('Oc1ccc(cc1)C(c2ccc(O)cc2)(C)C', 'ClC(=O)Cl'),
            ('ClC(=O)Oc1ccc(cc1)C(c2ccc(O)cc2)(C)C',),
        ),
        ## polycarbonate (non-phosgene route) - # BPA + diphenyl carbonate
        (
            "[#1]-[#8:1]-[!$([#6]=[#8]):2].*-[#8]-[#6:3](=[#8:4])-[*:5]>>[#8:1](-[!$([#6]=[#8]):2])-[#6:3](=[#8:4])-[*:5]",
            ('Oc1ccc(cc1)C(c2ccc(O)cc2)(C)C', 'O=C(Oc1ccccc1)Oc2ccccc2'),
            ('c1ccccc1OC(=O)Oc2ccc(cc2)C(C)(C)c3ccc(O)cc3',),
        ),
        ## polyurethane (isocyanate route) - # Bayer HDI + BDO
        (
            "[#1:1]-[#8:2]-[!$([#6]=[#8]):3].[#8:4]=[#6:5]=[#7:6]-[*:7]>>[#1:1]-[#7:6](-[#6:5](-[#8:2]-[!$([#6]=[#8]):3])=[#8:4])-[*:7]",
            ('O=C=N\CCCCCC/N=C=O', 'OCCCCO'),
            ('O=C=NCCCCCCNC(=O)OCCCCO',),
        ),
        ## polyurethane (non-isocyanate route) - # PCA (propylene carbonate acrylate) + hexamethylenediamine
        (
            "[*:1]-[#6:2]1(-[#1:3])-[#8:4]-[#6:5](=[#8:6])-[#8:7]-[#6:8]-1(-[*:9])-[#1:10].[#7:11](-[*:12])(-[#1:13])-[#1:14]>>[*:1]-[#6:2](-[#1:3])(-[#8:4]-[#1:13])-[#6:8](-[#8:7]-[#6:5](=[#8:6])-[#7:11](-[*:12])-[#1:14])(-[*:9])-[#1:10]",
            ('CC(=C)C(=O)OCC1COC(=O)O1', 'NCCCCCCN'),
            ('NCCCCCCNC(=O)OCC(-O)COC(=O)C(C)=C',),
        ),
        ## polyvinyl - polystyrene
        (
            "[*:1]-[#6:2](-[*:3])=[#6:4](-[#1:5])-[#1:6].[*:7]-[#6:8](=[#6:9](-[#1:10])-[#1:11])-[#1:12]>>[*:1]-[#6:2](-[*:3])(-[#6:4](-[#1:5])(-[#1:6])-[#6:8](-[*:7])=[#6:9](-[#1:10])-[#1:11])-[#1:12]",
            ('c1ccccc1C=C', 'c1ccccc1C=C'),
            ('c1ccccc1CCC(=C)c2ccccc2',),
        )
    ]
)
def test_reaction(rxn_smarts : Smarts, reactant_smiles : tuple[Smiles], expected_product_patterns : tuple[Smiles]) -> None:
    '''Test that a given reaction produces the correct products from known reactants'''
    rxn = AnnotatedReaction.from_smarts(rxn_smarts)
    reactants = rxn.valid_reactant_ordering(explicit_mols_from_SMILES(reactant_smiles))
    products_expected = explicit_mols_from_SMILES(expected_product_patterns)
    
    products_actual = rxn.react(reactants, keep_map_labels=False)
    assert len(products_actual) == len(products_expected), 'Mismatched number of products compared to expectation'

    for product_expected, product_actual in zip(products_expected, products_actual, strict=True):
        print(Chem.MolToSmiles(product_expected, canonical=True))
        print(Chem.MolToSmiles(product_actual, canonical=True))
        assert Chem.MolToSmiles(product_expected, canonical=True) == Chem.MolToSmiles(product_actual, canonical=True), 'Mismatches product'
    
    