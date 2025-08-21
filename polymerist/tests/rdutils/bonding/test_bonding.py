'''Behavioral tests for bond formation, dissolution, and splicing'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import pytest

from rdkit import Chem

from polymerist.rdutils.chemlabel import atom_idxs_by_map_numbers
from polymerist.rdutils.bonding.portlib import Port
from polymerist.rdutils.bonding.dissolution import dissolve_bond, decrease_bond_order
from polymerist.rdutils.bonding.formation import increase_bond_order
from polymerist.rdutils.bonding.substitution import saturate_ports
from polymerist.rdutils.bonding.permutation import swap_bonds


TEST_FLAVOR_PAIR : tuple[int, int] = (6, 9) # pair of linker favor labels to use throughout examples when they are called for
Port.bondable_flavors.insert(TEST_FLAVOR_PAIR) 

PATHOLOGICAL_FRAGMENTS : dict[str, str] = {
    'neutronium' : '[*:1]-[1*:2]',
    'ghost_water' : '[H:1]-[*:2]-[H:3]',
    'double_mid' : '[C:1](-[1*:2])(-[H:3])=[*:4]-[C:5](-[H:6])(-[H:7])(-[2*:8])'
}

@pytest.mark.parametrize(
    'smiles,map_num_pair,query_smarts',
    [
        (   
            '[1*:1]-[C:2](=[O:3])-[N:4](-[H:4])(-[H:5])',
            (2, 4), # split on single bond
            '[9*]-[NH2].[6*]C(=O)-[1*]',
        ),
        (   
            '[1*:1]-[C:2](=[O:3])-[N:4](-[H:4])(-[H:5])',
            (1, 2), # split on single bond to existing port, creating neutronium
            '[1*]-[6*].[9*:1]-C(=[O:3])-[NH2]',
        ),
        (   
            '[1*:1]-[C:2](=[O:3])-[N:4](-[H:4])(-[H:5])',
            (2, 3), # down-convert double bond to single bond
            '[1*:1]-C(-[6*])(-[O:3]-[9*])-[NH2]',
        ),
    ]
)
def test_decrease_bond_order(
        smiles : str,
        map_num_pair : tuple[int, int],
        query_smarts : str,
    ) -> None:
    '''Test that cutting bonds between atoms and replacing them with pairs of ports produces the expected molecular structure'''
    mol = Chem.RWMol(Chem.MolFromSmiles(smiles, sanitize=False))
    newmol = decrease_bond_order(mol, *atom_idxs_by_map_numbers(mol, *map_num_pair), new_flavor_pair=TEST_FLAVOR_PAIR)
    query = Chem.MolFromSmarts(query_smarts)

    assert newmol.HasSubstructMatch(query)

@pytest.mark.parametrize(
    'smiles,map_num_pair,query_smarts',
    [
        # NOTE: following 3 examples reverse the triple-bond dissolution in the test_decrease_bond_order() examples
        (
            '[N:1](-[6*:2])(-[6*:3])(-[6*:4]).[C:5](-[9*:6])(-[9*:7])(-[9*:8])[C:9](=[O:10])[O:11]-[H:12]',
            (1, 5),
            'N(-[6*])(-[6*])C(-[9*])(-[9*])C(=O)[OH]',
        ),
        (
            '[N:1]([6*:2])([6*:4])[C:5]([9*:6])([9*:7])[C:9](=[O:10])[O:11][H:12]',
            (1, 5),
            'N(-[6*])=C(-[9*])C(=O)[OH]',
        ),
        (
            '[N:1]([6*:4])=[C:5]([9*:7])[C:9](=[O:10])[O:11][H:12]',
            (1, 5),
            '[N:1]#[C:2]-[C:3](=[O:4])-[O:5]-[H:6]',
        )
    ]
)
def test_increase_bond_order(
        smiles : str,
        map_num_pair : tuple[int, int],
        query_smarts : str,
    ):
    '''Test that forming new bonds from pairs of compatible ports produces the expected molecular structure'''
    mol = Chem.RWMol(Chem.MolFromSmiles(smiles, sanitize=False))
    newmol = increase_bond_order(mol, *atom_idxs_by_map_numbers(mol, *map_num_pair))
    query = Chem.MolFromSmarts(query_smarts)
    
    assert newmol.HasSubstructMatch(query)
    
@pytest.mark.parametrize(
    'smiles,map_num_pair,query_smarts',
    [
        (
            '[1*:1]=[C:2](-[H:6])-[N:3](-[H:4])-[H:5]',
            (2, 3),
            '[9*]-[NH2].[6*]-C(-[H])=[*]',
        ),
        # NOTE: below examples are of a made-up molecule cooked up to contain at least one single, double, and triple bond
        (
            '[N:1]#[C:2]-[C:3](=[O:4])-[O:5]-[H:6]',
            (2, 3), # dissolve single bond
            'N#C-[6*].[9*]-C(=O)[OH]',
        ),                                
        (
            '[N:1]#[C:2]-[C:3](=[O:4])-[O:5]-[H:6]',
            (3, 4), # dissolve double bond
            'N#CC(-[6*])(-[6*])[OH].[9*]-O-[9*]',
        ),                     
        (
            '[N:1]#[C:2]-[C:3](=[O:4])-[O:5]-[H:6]',
            (1, 2), # dissolve triple bond
            'N(-[6*])(-[6*])(-[6*]).C(-[9*])(-[9*])(-[9*])C(=O)[OH]',
        ), 
    ]
)
def test_dissolve_bond(
        smiles : str,
        map_num_pair : tuple[int, int],
        query_smarts : str,
    ) -> None:
    '''Test that dissolving bonds completely into component Ports produces the expected molecular structure'''
    mol = Chem.RWMol(Chem.MolFromSmiles(smiles, sanitize=False))
    newmol = dissolve_bond(mol, *atom_idxs_by_map_numbers(mol, *map_num_pair), new_flavor_pair=TEST_FLAVOR_PAIR)
    query = Chem.MolFromSmarts(query_smarts)
    
    assert newmol.HasSubstructMatch(query)
          
@pytest.mark.parametrize(
    'core_smiles,cap_smiles,query_smarts',
    [
        ( # splice single bonds
            'C(-[6*])(-[6*])(-[6*])(-[6*])',
            '[9*]-[O]-[H]',
            'C([OH])([OH])([OH])([OH])'
        ),
        ( # splice single bonds - selective
            'C(-[6*])(-[6*])(-[5*])(-[5*])',
            '[9*]-[O]-[H]',
            'C([OH])([OH])(-[5*])(-[5*])'
        ),
        ( # splice double bonds
            '[6*]=C=[6*]',
            '[9*]=O',
            'O=C=O',
        ),
        ( # splice double bonds - selective
            '[6*]=C=[4*]',
            '[9*]=O',
            'O=C=[4*]',
        ),
        ( # splice triple bonds
            '[6*]#CC#[6*]',
            '[9*]#N',
            'N#CC#N',
        ),
        ( # splice triple bonds - SELECTIVE
            '[6*]#CC#[5*]',
            '[9*]#N',
            'N#CC#[*5]',
        ),
    ]
) # DEV: implicitly also tests substitution.splice_atoms()
def test_saturate_ports(
        core_smiles : str,
        cap_smiles : str,
        query_smarts : str,
    ):
    '''Test that capping all available ports on a "core" molecule with terminal "cap" molecules yields the expected structure'''
    core_mol = Chem.MolFromSmiles(core_smiles, sanitize=False)
    cap_mol = Chem.MolFromSmiles(cap_smiles, sanitize=False)
    query = Chem.MolFromSmarts(query_smarts)

    newmol = saturate_ports(core_mol, cap=cap_mol, flavor_to_saturate=TEST_FLAVOR_PAIR[0])
    assert newmol.HasSubstructMatch(query)
    
@pytest.mark.parametrize(
    'smiles,derangement,query_smarts',
    [
        ( # polyamide condensation
            '[*:1]-[C:2](=[O:3])-[O:4]-[H:5].[*:6]-[N:7](-[H:8])-[H:9]',
            {7 : (9, 2), 4 : (2, 9)},
            '[H]-O-[H].*-N(-[H])C(=O)-*',
        ),
        ( # polycarbonate (phosgene route) rxn definition
            '[*:1]-[O:2]-[H:3].[*:4]-[C:5](=[O:6])-[Cl:7]',
            {5 : (7, 2), 3 : (2, 7)},
            '[H]-[Cl].*-C(=O)O-*',
        ),
        ( # polyurethane (isocyanate) rxn definition
            '[*:1]-[O:2]-[H:3].[*:4]-[N:5]=[C:6]=[O:7]',
            {2 : (3, 6), 5 : (6, 3)},
            '*-O-C(=O)-N(-[H])-*',
        ),
    ]
)
def test_swap_bonds(
        smiles : str,
        derangement : dict[int, tuple[int, int]], 
        query_smarts : str
    ) -> None:
    '''Test that permuting bonds within a Mol yields the expected structure'''
    mol = Chem.RWMol(Chem.MolFromSmiles(smiles, sanitize=False))
    newmol = swap_bonds(mol, derangement, show_steps=False, in_place=False)
    query = Chem.MolFromSmarts(query_smarts)
    
    assert newmol.HasSubstructMatch(query)
    