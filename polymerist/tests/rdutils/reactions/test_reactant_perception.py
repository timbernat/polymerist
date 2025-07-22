'''Behavioral tests for reactant chemical perception for `reactions.reactions` submodule'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import pytest

from rdkit import Chem
from polymerist.rdutils.reactions.reactions import AnnotatedReaction
from . import RXNS


# DEFINING REACTANTS
@pytest.fixture(scope='module')
def test_reactant_pool_smiles() -> list[str]:
    return [
        'CC(O)C(Cl)C(Cl)CC(=O)O',
        'c1ccc(Cl)cc1O',
        'c1ccc(C(=O)Cl)cc1C(=O)O',
        'FC(F)(F)F',
        'c1(O)c(O)c(C(=O)Cl)c(O)c(C(=O)Cl)c1',
        'O=C=N\\CCCCCC/N=C=O',
        'OCCC(N=C=O)CCO',
        'C=CC(=O)OCC1COC(=O)O1',
        'NCCCCCCN',
        'NCCC(C(=O)O)CCCN',
        'C=Cc1ccccc1',
        'C=C(Cl)',
        'FC(F)=C',
        'O1C(=O)C(=O)C(N)(N)C(=O)C(=O)1',
        # NOTE: have deliberately excluded compatible reactants for polyimides to test whether negative results are correctly returned during reactant perception
    ]
    
@pytest.fixture(scope='module')
def test_reactant_pool(test_reactant_pool_smiles : list[str]) -> list[Chem.Mol]:
    TEST_REACTANTS : list[Chem.Mol] = []
    for smiles in test_reactant_pool_smiles:
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
        mol = Chem.AddHs(mol, addCoords=True)
        Chem.Kekulize(mol, clearAromaticFlags=True)
        Chem.SetAromaticity(mol, model=Chem.AROMATICITY_MDL)
        Chem.SanitizeMol(mol)
        
        TEST_REACTANTS.append(mol)
    return TEST_REACTANTS

FUNCTIONAL_GROUP_IDXS : dict[str, set[int]] = { # indices of which test mols contain each type of functional group
    'amine' : {8, 9, 13},
    'acyl_chloride' : {2, 4},
    'hydroxyl' : {0, 1, 4, 6},
    'carboxyl' : {0, 2, 9},
    'cyclocarbonate' : {7},
    'isocyanate' : {5, 6}, 
    'vinyl' : {7, 10, 11},
    'anhydride' : {13},
}

    
# TESTING REACTANT PERCEPTION
@pytest.mark.parametrize('rxn', RXNS.values())
def test_reactant_perception_affirmative(rxn : AnnotatedReaction, test_reactant_pool : list[Chem.Mol]) -> None:
    '''Test that a reaction correctly detects when a pool of reactants contains at least one reactable subset (without enumerating those reactants)'''
    assert rxn.has_reactable_subset(test_reactant_pool, allow_resampling=False) # test reactants above are defined to participate in at least one of the example rxn mechanisms

@pytest.mark.parametrize(
    'rxn,idxs_to_exclude',
    [
        (RXNS['polyvinyl'], FUNCTIONAL_GROUP_IDXS['vinyl']),
        (RXNS['polycarbonate'], FUNCTIONAL_GROUP_IDXS['acyl_chloride']),
        (RXNS['polyimide'], FUNCTIONAL_GROUP_IDXS['anhydride']),
        (RXNS['polyamide'], FUNCTIONAL_GROUP_IDXS['amine']),
        (RXNS['polyester'], FUNCTIONAL_GROUP_IDXS['carboxyl']),
        (RXNS['polyurethane_isocyanate'], FUNCTIONAL_GROUP_IDXS['isocyanate']),
        (RXNS['polyurethane_nonisocyanate'], FUNCTIONAL_GROUP_IDXS['amine']),
    ]
)
def test_reactant_perception_negative(
        rxn : AnnotatedReaction, 
        test_reactant_pool : list[Chem.Mol], 
        idxs_to_exclude : set[int], 
    ) -> None:
    '''Test that a reaction correctly detects when a pool of reactants DOES NOT contain at least one reactable subset (without enumerating those reactants)'''
    reactants_pool = [ # simplifies checking if removing a certain subset of reactants changes whether a mechanism should be possible
        reactant
            for i, reactant in enumerate(test_reactant_pool)
                if (i not in idxs_to_exclude)
    ]
    assert not rxn.has_reactable_subset(reactants_pool, allow_resampling=False) # by construction, each of these reactions should be chemically impossible

@pytest.mark.parametrize(
    'rxn,reactive_indices_expected',
    [ # all ordered pairs of reactants compatible with each reaction
        (RXNS['polyamide'], {(8, 0), (9, 0), (8, 2), (9, 2), (8, 9), (13, 0), (13, 2), (13, 9)}),
        (RXNS['polyester'], {(0, 2), (0, 9), (1, 0), (1, 2), (1, 9), (4, 0), (4, 2), (4, 9), (6, 0), (6, 2), (6, 9)}),
        (RXNS['polyvinyl'], {(7, 11), (7, 10), (10, 11), (10, 7), (11, 10), (11, 7), (12, 11), (12, 10), (12, 7)}),
        (RXNS['polycarbonate'], {(0, 2), (0, 4), (1, 2), (1, 4), (4, 2), (6, 2), (6, 4)}),
        (RXNS['polyimide'], {(8, 13), (9, 13)}),
        (RXNS['polyurethane_isocyanate'], {(0, 5), (0, 6), (1, 5), (1, 6), (4, 5), (4, 6), (6, 5)}),
        (RXNS['polyurethane_nonisocyanate'], {(7, 8), (7, 9), (7, 13)}),
    ]
)
def test_reactant_order_enumeration_separated(rxn : AnnotatedReaction, test_reactant_pool : list[Chem.Mol], reactive_indices_expected : set[tuple[int, int]]) -> None:
    '''
    Test that reactions are able to enumerate all compatible reactants from a collection of arbitrary reactants,
    where each reactant is only allowed to be used for at most 1 functional group'''
    assert set(
        reactant_ordering 
            for reactant_ordering in rxn.enumerate_valid_reactant_orderings(test_reactant_pool, allow_resampling=False, as_mols=False)
    ) == reactive_indices_expected
        
@pytest.mark.parametrize(
    'rxn,autoreactant_indices_expected',
    [ # indices of all molecules which can conduct a reaction by themselves (i.e. pure-component)
        (RXNS['polyamide'], {9}),
        (RXNS['polyester'], {0}),
        (RXNS['polyvinyl'], {7, 10, 11}),
        (RXNS['polycarbonate'], {4}),
        (RXNS['polyimide'], {13}),
        (RXNS['polyurethane_isocyanate'], {6}),
        (RXNS['polyurethane_nonisocyanate'], set()),
    ]
)
def test_reactant_order_enumeration_autoreactions(rxn : AnnotatedReaction, test_reactant_pool : list[Chem.Mol], autoreactant_indices_expected : set[int]) -> None:
    '''
    Test that intramolecular reactions are identified (test of resampling ability, with single-molecule input sets)'''
    assert {
        reactant_idx
            for reactant_idx, reactant in enumerate(test_reactant_pool)
                if rxn.valid_reactant_ordering([reactant], as_mols=False, allow_resampling=True) is not None
    } == autoreactant_indices_expected