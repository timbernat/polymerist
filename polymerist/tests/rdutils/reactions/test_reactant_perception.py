'''Behavioral tests for reactant chemical perception for `reactions.reactions` submodule'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import pytest

from rdkit import Chem
from polymerist.rdutils.reactions.reactions import AnnotatedReaction


# DEFINING REACTIONS
RXN_SMARTS : dict[str, str] = {
    'polyamide'     : '[#7:1](-[*:2])(-[#1])-[#1:3].[#8](-[#6:4](=[#8:5])-[*:6])-[#1]>>[#7:1](-[*:2])(-[#1:3])-[#6:4](=[#8:5])-[*:6]',
    'polyvinyl'     : '[*:1]-[#6:2](-[*:3])=[#6:4](-[#1:5])-[#1:6].[*:7]-[#6:8](=[#6:9](-[#1:10])-[#1:11])-[#1:12]>>[*:1]-[#6:2](-[*:3])(-[#6:4](-[#1:5])(-[#1:6])-[#6:8](-[*:7])=[#6:9](-[#1:10])-[#1:11])-[#1:12]',
    'polycarbonate' : '[#1]-[#8:1]-[!$([#6]=[#8]):2].[#17]-[#6:3](=[#8:4])-[*:5]>>[#8:1](-[!$([#6]=[#8]):2])-[#6:3](=[#8:4])-[*:5]',
    'polyurethane_isocyanate': '[#1:1]-[#8:2]-[!$([#6]=[#8]):3].[#8:4]=[#6:5]=[#7:6]-[*:7]>>[#1:1]-[#7:6](-[#6:5](-[#8:2]-[!$([#6]=[#8]):3])=[#8:4])-[*:7]',
    # NOTE: have deliberately excluded compatible reactants for polyimides to test whether negative results are correctly returned during reactant perception
    'polyimide' : '[#7:1](-[*:2])(-[#1])-[#1].[*:3]-[#6:4](=[#8:5])-[#8]-[#6:6](=[#8:7])-[*:8]>>[#7:1](-[*:2])(-[#6:4](-[*:3])=[#8:5])-[#6:6](=[#8:7])-[*:8]',
}
RXNS : dict[str, AnnotatedReaction] = { # TODO: write tests for initialization (i.e. from_smarts() and from_rdmols)
    rxn_name : AnnotatedReaction.from_smarts(rxn_smarts)
        for rxn_name, rxn_smarts in RXN_SMARTS.items()
}

# DEFINING REACTANTS
@pytest.fixture(scope='module')
def test_reactant_smiles() -> list[str]:
    return [ # for polycarbonate
        'CC(O)C(Cl)C(Cl)CC(=O)O',
        'c1ccc(Cl)cc1O',
        'c1ccc(C(=O)Cl)cc1C(=O)O',
        'FC(F)(F)F',
        'c1(O)c(O)c(C(=O)Cl)c(O)c(C(=O)Cl)c1',
        "O=C=N\\CCCCCC/N=C=O",
        "OCCC(N=C=O)CCO",
        'C=CC(=O)OCC1COC(=O)O1',
        'NCCCCCCN',
        'NCCC(C(=O)O)CCCN',
        'C=Cc1ccccc1',
        'C=C(Cl)',
        'FC(F)=C',
    ]
    
@pytest.fixture(scope='module')
def test_reactants(test_reactant_smiles : list[str]) -> list[Chem.Mol]:
    TEST_REACTANTS : list[Chem.Mol] = []
    for smiles in test_reactant_smiles:
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
        mol = Chem.AddHs(mol, addCoords=True)
        Chem.Kekulize(mol, clearAromaticFlags=True)
        Chem.SetAromaticity(mol, model=Chem.AROMATICITY_MDL)
        Chem.SanitizeMol(mol)
        
        TEST_REACTANTS.append(mol)
    return TEST_REACTANTS


# TESTING REACTANT PERCEPTION
reactant_pair_indices : dict[AnnotatedReaction, set[tuple[int]]] = { # all ordered pairs of reactants compatible with each reaction
    RXNS['polyamide']               : {(8, 0), (9, 0), (8, 2), (9, 2), (8, 9)},
    RXNS['polyvinyl']               : {(7, 11), (7, 10), (10, 11), (10, 7), (11, 10), (11, 7), (12, 11), (12, 10), (12, 7)},
    RXNS['polycarbonate']           : {(0, 2), (0, 4), (1, 2), (1, 4), (4, 2), (6, 2), (6, 4)},
    RXNS['polyurethane_isocyanate'] : {(0, 5), (0, 6), (1, 5), (1, 6), (4, 5), (4, 6), (6, 5)},
}
autoreactant_indices : dict[AnnotatedReaction, set[tuple[int]]] = { # indices of all molecules which can conduct a reaction by themselves (i.e. pure-component)
    RXNS['polyamide']               : {9},
    RXNS['polyvinyl']               : {7, 10, 11},
    RXNS['polycarbonate']           : {4},
    RXNS['polyurethane_isocyanate'] : {6},
}

@pytest.mark.parametrize(
    'rxn,has_reactable_subset_expected',
    [(RXNS['polyamide'], True), (RXNS['polyvinyl'], True), (RXNS['polycarbonate'], True), (RXNS['polyurethane_isocyanate'], True), (RXNS['polyimide'], False)]
)
def test_reactant_perception(rxn : AnnotatedReaction, test_reactants : list[Chem.Mol], has_reactable_subset_expected : bool) -> None:
    '''Test that a reaction can detect if a set of reactants contains at least one reactable subset (not yet enumerating those reactants)'''
    assert rxn.has_reactable_subset(test_reactants) == has_reactable_subset_expected

@pytest.mark.parametrize('rxn,reactive_indices_expected', reactant_pair_indices.items())
def test_reactant_order_enumeration_separated(rxn : AnnotatedReaction, test_reactants : list[Chem.Mol], reactive_indices_expected : set[tuple[int, int]]) -> None:
    '''
    Test that reactions are able to enumerate all compatible reactants from a collection of arbitrary reactants,
    where each reactant is only allowed to be used for at most 1 functional group'''
    assert set(
        reactant_ordering 
            for reactant_ordering in rxn.enumerate_valid_reactant_orderings(test_reactants, allow_resampling=False, as_mols=False)
    ) == reactive_indices_expected
        
@pytest.mark.parametrize('rxn,autoreactant_indices_expected', autoreactant_indices.items())
def test_reactant_order_enumeration_autoreactions(rxn : AnnotatedReaction, test_reactants : list[Chem.Mol], autoreactant_indices_expected : set[int]) -> None:
    '''
    Test that intramolecular reactions are identified (test of resampling ability, with single-molecule input sets)'''
    assert {
        i
            for i, reactant in enumerate(test_reactants)
                if rxn.valid_reactant_ordering([reactant], as_mols=False, allow_resampling=True) is not None
    } == autoreactant_indices_expected