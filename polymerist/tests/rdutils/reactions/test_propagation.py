'''Behavioral tests to ensure reaction fragment enumeration behaves as expected'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import pytest

from rdkit import Chem
from rdkit.Chem.rdmolops import AromaticityModel, SANITIZE_ALL, AROMATICITY_MDL

from . import RXNS

from polymerist.smileslib.cleanup import expanded_SMILES
from polymerist.rdutils.sanitization import sanitize_mol, explicit_mol_from_SMILES

from polymerist.rdutils.reactions.reactions import AnnotatedReaction
from polymerist.rdutils.reactions.reactors import PolymerizationReactor
from polymerist.rdutils.reactions.fragment import (
    IBIS,
    ReseparateRGroups,
    CutMinimumCostBondsStrategy,
)

# CONFIG CONSTANTS
@pytest.fixture(scope='module')
def fragment_strategy() -> IBIS:
    '''Fragmentation strategy to use for the tests'''
    return CutMinimumCostBondsStrategy()

@pytest.fixture(scope='module')
def aromaticity_model() -> IBIS:
    '''Aromaticity model to use when sanitizing Mols'''
    # NOTE: polyimide tests WILL fail if aromaticity is not MDL (e.g. is "AROMATICITY_RDKIT")
    return AROMATICITY_MDL 
    
    
# DEFINING REACTANTS  
@pytest.mark.parametrize(
    'reactant_smiles,rxn,allow_resampling,rxn_depth_max, fragment_smarts_expected',
    (
        ## common condensation reactions - one-step only
        ( ### polyamide
            {
                'NCCCCCC(=O)O',
                'NCCCCCCN',
                'O=C(O)c1ccc(C(=O)O)cc1',    
            },
            RXNS['polyamide'], False, 1, 
            {
                'NCCCCCCN*',
                '*NCCCCCC(=O)O',
                'NCCCCCC(=O)*',
                'O=C(*)c1ccc(C(=O)O)cc1',
                'NCCCCCC(=O)O',
                'NCCCCCCN',
                'O=C(O)c1ccc(C(=O)O)cc1',      
            },
        ),
        ( ### polycarbonate
            {
                'CC(C)(c1ccc(O)cc1)c1ccc(O)cc1',
                'O=C(Cl)Cl',
                'OCCCCO',
            },
            RXNS['polycarbonate'], False, 1,
            {
                'O=C(Cl)Cl',
                '*C(=O)Cl',
                'OCCCCO',
                '*OCCCCO',
                'CC(C)(c1ccc(O)cc1)c1ccc(O)cc1',
                '*Oc1ccc(C(C)(C)c2ccc(O)cc2)cc1',
            },
        ),
        ( ### polyimide - PMDA + ODA (Kapton), for testing multi-ring cuts and aromaticity setting
            {
                'c1cc(N)ccc1Oc1ccc(N)cc1',
                'c1c2C(=O)OC(=O)c2cc3C(=O)OC(=O)c31',
            },
            RXNS['polyimide'], False, 1,
            {
                'c1cc(N)ccc1Oc1ccc(N)cc1',
                'c1cc(*)ccc1Oc1ccc(N)cc1',
                'c1c2C(=O)OC(=O)c2cc3C(=O)OC(=O)c31',
                'c1c2C(=O)N(*)C(=O)c2cc3C(=O)OC(=O)c31',
            },
        ),
        ( ### polyester - implicitly also tests that hydroxyls and carboxyls are correctly differentiated in the rxn definition
            {
                'O=C(O)CCC(=O)O',
                'OCCCCO',
                'Oc1ccc(C(=O)O)o1',
            },
            RXNS['polyester'], False, 1,
            {
                'OCCCCO',
                '*OCCCCO',
                'O=C(O)CCC(=O)O',
                '*C(=O)CCC(=O)O',
                'O=C(O)c1ccc(O)o1',
                '*Oc1ccc(C(=O)O)o1',
                '*C(=O)c1ccc(O)o1',
            },
        ),
        ## common condensation reactions - full rxn tree search
        ( ### polyamide
            {
                'NCCCCCC(=O)O',
                'NCCCCCCN',
                'O=C(O)c1ccc(C(=O)O)cc1',    
            },
            RXNS['polyamide'], False, 5, 
            {
                '*C(=O)CCCCCN',
                '*C(=O)c1ccc(C(*)=O)cc1',
                '*C(=O)c1ccc(C(=O)O)cc1',
                '*NCCCCCC(*)=O',
                '*NCCCCCC(=O)O',
                '*NCCCCCCN',
                '*NCCCCCCN*',
                'NCCCCCC(=O)O',
                'NCCCCCCN',
                'O=C(O)c1ccc(C(=O)O)cc1',  
            },
        ),
        ( ### polycarbonate
            {
                'CC(C)(c1ccc(O)cc1)c1ccc(O)cc1',
                'O=C(Cl)Cl',
                'OCCCCO',
            },
            RXNS['polycarbonate'], False, 5,
            {
                'O=C(Cl)Cl',
                '*C(=O)Cl',
                '*C(*)=O',
                'OCCCCO',
                '*OCCCCO',
                '*OCCCCO*',
                'CC(C)(c1ccc(O)cc1)c1ccc(O)cc1',
                '*Oc1ccc(C(C)(C)c2ccc(O)cc2)cc1',
                '*Oc1ccc(C(C)(C)c2ccc(O*)cc2)cc1',
            },
        ),
        ( ### polyimide - PMDA + ODA (Kapton), for testing multi-ring cuts and aromaticity setting
            {
                'c1cc(N)ccc1Oc1ccc(N)cc1',
                'c1c2C(=O)OC(=O)c2cc3C(=O)OC(=O)c31',
            },
            RXNS['polyimide'], False, 5,
            {
                'c1cc(N)ccc1Oc1ccc(N)cc1',
                'c1cc(*)ccc1Oc1ccc(N)cc1',
                'c1cc(*)ccc1Oc1ccc(*)cc1',
                'c1c2C(=O)OC(=O)c2cc3C(=O)OC(=O)c31',
                'c1c2C(=O)N(*)C(=O)c2cc3C(=O)OC(=O)c31',
                'c1c2C(=O)N(*)C(=O)c2cc3C(=O)N(*)C(=O)c31',
            },
        ),
        ( ### polyester - implicitly also tests that hydroxyls and carboxyls are correctly differentiated in the rxn definition
            {
                'O=C(O)CCC(=O)O',
                'OCCCCO',
                'Oc1ccc(C(=O)O)o1',
            },
            RXNS['polyester'], False, 5,
            {
                'OCCCCO',
                '*OCCCCO',
                '*OCCCCO*',
                'O=C(O)CCC(=O)O',
                '*C(=O)CCC(=O)O',
                '*C(=O)CCC(*)=O',
                'O=C(O)c1ccc(O)o1',
                '*Oc1ccc(C(=O)O)o1',
                '*C(=O)c1ccc(O)o1',
                '*Oc1ccc(C(*)=O)o1',
            },
        ),
        ## autopolymerization, WITHOUT resampling - should halt immediately
        ( ### polyamide - 6-aminohexanoic acid
            {'NCCCCCC(=O)O'},
            RXNS['polyamide'], False, 3,
            {'NCCCCCC(=O)O'},
        ),
        ( ### vinyl - polystyrene
            {'C=Cc1ccccc1'},
            RXNS['polyvinyl'], False, 3,
            {'C=Cc1ccccc1'},
        ),
        ## autopolymerization, WITH resampling - should produce new fragments
        ( ### polyamide - 6-aminohexanoic acid
            {'NCCCCCC(=O)O'},
            RXNS['polyamide'], True, 3,
            {
                '*NCCCCCC(*)=O',
                '*C(=O)CCCCCN',
                '*NCCCCCC(=O)O',
                'NCCCCCC(=O)O',
            }
        ),
        ( ### vinyl - polystyrene
            {'C=Cc1ccccc1'}, 
            RXNS['polyvinyl'], True, 3,
            {
            '*CC(*)c1ccccc1',
            '*CCc1ccccc1',
            '*C(=C)c1ccccc1',
            'C=Cc1ccccc1',
        }),
    )
)

# TESTS PROPER
def test_propagation_pooled(
    reactant_smiles : set[str],
    rxn : AnnotatedReaction,
    allow_resampling : bool,
    rxn_depth_max : int,
    fragment_smarts_expected : set[str],
    fragment_strategy : IBIS,
    aromaticity_model : AromaticityModel,
) -> None:
    '''Test that propagation of monomers produces the expected fragments'''
    reactants : list[Chem.Mol] = [
        explicit_mol_from_SMILES(
            smiles,
            assign_map_nums=False,
            sanitize_ops=SANITIZE_ALL,
            # NOTE: polyimide examples WILL fail if aromaticity is not MDL (e.g. "AROMATICITY_RDKIT")!
            aromaticity_model=aromaticity_model,
        )
            for smiles in reactant_smiles
    ]
    fragment_smarts_canon : set[str] = set(
        Chem.CanonSmiles(smiles_expected)
            for smiles_expected in fragment_smarts_expected
    )
        
    reactor = PolymerizationReactor(rxn, fragment_strategy=fragment_strategy)
    fragments : dict[str, Chem.Mol] = reactor.propagate_pooled(
        reactants,
        rxn_depth_max=rxn_depth_max,
        allow_resampling=allow_resampling,
        aromaticity_model=aromaticity_model,
    )
    
    assert (set(fragments.keys()) == fragment_smarts_canon)