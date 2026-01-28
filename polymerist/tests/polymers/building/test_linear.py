'''Tests construction of structures for linear copolymers (and relevant subfamilies, e.g. homopolymers)'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import pytest

from collections import Counter
from itertools import combinations, product as cartesian_product

from polymerist.rdutils.sanitization import explicit_mol_from_SMILES

from polymerist.polymers.monomers.repr import MonomerGroup
from polymerist.polymers.monomers.fragments import PE_FRAGMENTS, MPD_TMC_FRAGMENTS, PEG_PLGA_FRAGMENTS
from polymerist.polymers.building.linear import build_linear_polymer
from polymerist.polymers.building.mbconvert import mbmol_to_rdmol
from polymerist.polymers.exceptions import MorphologyError, PartialBlockSequence, EmptyBlockSequence


# FIXTURES AND INPUT SETUP FOR TESTS
@pytest.fixture(scope='function')
def monogrp_polyethylene() ->  MonomerGroup:
    return MonomerGroup(monomers=PE_FRAGMENTS)

@pytest.fixture(scope='function')
def monogrp_mpd_tmc() ->  MonomerGroup:
    return MonomerGroup(monomers=MPD_TMC_FRAGMENTS)

@pytest.fixture(scope='function')
def monogrp_peg_plga() ->  MonomerGroup:
    return MonomerGroup(monomers=PEG_PLGA_FRAGMENTS)

@pytest.fixture(scope='function')
def monogrp_halogen_marked() -> MonomerGroup:
    '''
    A MonomerGroup with single-carbon repeat units marked by halogens
    to make repeat units easy to distinguish when testing repeat unit sequencing
    '''
    monogrp = MonomerGroup()
    mono_smarts = {
        # fluorines
        'fluor_term_1' : 'OC(F)-*',
        'fluor_mid_1'  : '*-C(F)-*',
        'fluor_mid_2'  : '*-C(F)(F)-*',
        'fluor_term_2' : '*-C(F)C(=O)O',
        # chlorines
        'chlor_term_1' : 'OC(Cl)-*',
        'chlor_mid_1'  : '*-C(Cl)-*',
        'chlor_mid_2'  : '*-C(Cl)(Cl)-*',
        'chlor_term_2' : '*-C(Cl)C(=O)O',
        # bromines
        'brom_term_1' : 'OC(Br)-*',
        'brom_mid_1'  : '*-C(Br)-*',
        'brom_mid_2'  : '*-C(Br)(Br)-*',
        'brom_term_2' : '*-C(Br)C(=O)O',
    }
    # term_orient deliberately omitted - to be set in test
    for resname, smarts in mono_smarts.items():
        monogrp.add_monomer(resname, smarts)

    return monogrp

SEQS_TO_HALOGEN_KERNEL : dict[tuple[str, ...], str] = {
    ('A', 'S', '$') : 'C(F)' , # for homopolymer, ANY single symbol should given same behavior
    ('AA', 'BB', '$$') : 'C(F)C(F)' , # BB is confusing, but point is this is still a homopolymer, so symbol choice still doesn't matter
    ('BAB', 'qaq', '101', '212') : 'C(F)(F)C(F)C(F)(F)',
    ('BADC', 'badc', '2143') : 'C(F)(F)C(F)C(Cl)(Cl)C(Cl)',
    ('BACC', 'bacc', '1022') : 'C(F)(F)C(F)C(Cl)C(Cl)',
    ('ABCFED', 'acezyx', '012543', '013964') : 'C(F)C(F)(F)C(Cl)C(Br)(Br)C(Br)C(Cl)(Cl)',
}

TERM_ORIENT_TO_SMILES_HALOGEN : dict[str, dict[str, str]] = {
    # fluorines
    'fluor_term_1' : {
        'head' : 'OC(F)',
        'tail' : 'C(F)O',
    },
    'fluor_term_2' : {
        'head' : 'OC(=O)C(F)',
        'tail' : 'C(F)C(=O)O',
    },
    # chlorines
    'chlor_term_1' : {
        'head' : 'OC(Cl)',
        'tail' : 'C(Cl)O',
    },
    'chlor_term_2' : {
        'head' : 'OC(=O)C(Cl)',
        'tail' : 'C(Cl)C(=O)O',
    },
    # bromines
    'brom_term_1' : {
        'head' : 'OC(Br)',
        'tail' : 'C(Br)O',
    },
    'brom_term_2' : {
        'head' : 'OC(=O)C(Br)',
        'tail' : 'C(Br)C(=O)O',
    },
}

def compile_sequence_test_inputs(kernel_repeats : tuple[int, ...]) -> list[tuple[str, dict[str, str], int, str]]:
    '''Helper methods for compiling flattened inputs for linear polymer builder sequence testing'''
    # pick out term group orientation dicts and accompaying SMILES caps
    term_orients_and_smiles : list[tuple[dict[str, str], tuple[str, str]]] = [
        ( # special case for unspecified term orients - inferred to be first two terminal groups
            dict(),
            (
                TERM_ORIENT_TO_SMILES_HALOGEN['fluor_term_1']['head'],
                TERM_ORIENT_TO_SMILES_HALOGEN['fluor_term_2']['tail'],
            ),
        )
    ]

    TERM_ORIENT_PAIRS : tuple[tuple[str, str], ...] = (
        ('fluor_term_1', 'fluor_term_2'),
        ('chlor_term_1', 'chlor_term_2'),
        ('chlor_term_2', 'chlor_term_1'),
        ('fluor_term_2', 'brom_term_1'),
        ('chlor_term_1', 'brom_term_1'),
    )
    # for (runame_head, runame_tail) in combinations(TERM_ORIENT_TO_SMILES_HALOGEN, 2): # opting to remove ALL combinations to avoid ballooning number of tests
    for (runame_head, runame_tail) in TERM_ORIENT_PAIRS: # opting to remove ALL combinations to avoid ballooning number of tests
        term_orients : dict[str, str] = {
            'head' : runame_head,
            'tail' : runame_tail,
        }
        term_smiles : tuple[str, str] = (
            TERM_ORIENT_TO_SMILES_HALOGEN[runame_head]['head'],
            TERM_ORIENT_TO_SMILES_HALOGEN[runame_tail]['tail'],
        )
        term_orients_and_smiles.append( (term_orients, term_smiles) )

    # pick out individual test inputs from cartesian product over options
    sequence_test_inputs : list[tuple[str, dict[str, str], int, str]] = []
    for ((sequences, kernel), (term_orient, term_smiles), n_kernel_repeats) in cartesian_product(
        SEQS_TO_HALOGEN_KERNEL.items(),
        term_orients_and_smiles,
        kernel_repeats,
    ):
        head_smiles, tail_smiles = term_smiles
        smiles_expected = ''.join([head_smiles] + n_kernel_repeats*[kernel] + [tail_smiles])

        for sequence in sequences:
            n_monomers = 2 + len(sequence)*n_kernel_repeats # 2 accounts for end groups
            sequence_test_inputs.append( (sequence, term_orient, n_monomers, smiles_expected) )

    return sequence_test_inputs


# TESTS PROPER
@pytest.mark.parametrize(
    'monomers, term_orient, n_monomers, sequence, minimize_sequence, allow_partial_sequences, energy_minimize',
    [
        # Polyethylene
        ('monogrp_polyethylene', {}, 7, 'A', True, True, False), # test end group autogen (should only have 1 term group)
        ('monogrp_polyethylene', {'head':'PE1', 'tail' : 'PE1'}, 7, 'A', True, True, False), # test explicit head-tail (result here should be different from autogen structure)
        ('monogrp_polyethylene', {'head':'PE1', 'tail' : 'PE1'}, 7, 'A', True, False, False), # test partial sequences (irrelevant here)
        ('monogrp_polyethylene', {'head':'PE1', 'tail' : 'PE1'}, 7, 'A', True, False, True), # test energy minimization doesn't crash
        pytest.param(  # will fail due to too few monomers for given sequence - 
            'monogrp_polyethylene', {}, 7, '', True, True, False, # NOTE: need to have partials enabled, since failure happens ONLY once sequence is passed to mbuild
            marks=pytest.mark.xfail(
                raises=EmptyBlockSequence,
                reason='Sequence provided must be nonempty',
                strict=True,
            )
        ),
        pytest.param(  # will fail due to too few monomers for given sequence - 
            'monogrp_polyethylene', {'head':'PE1', 'tail' : 'PE1'}, 7, 'AB', True, True, False, # NOTE: need to have partials enabled, since failure happens ONLY once sequence is passed to mbuild
            marks=pytest.mark.xfail(
                raises=ValueError,
                reason='Fewer unique monomers defined than called for by target sequence',
                strict=True,
            )
        ),
        # MPD-TMC
        ('monogrp_mpd_tmc', {'head':'MPD-1', 'tail' : 'TMC-1'}, 8, 'A', True, True, False), # correctly-specified: explicit end groups, only linear middle monomers, and whole number of sequence repeats
        pytest.param(
            'monogrp_mpd_tmc', {'head':'MPD-1', 'tail' : 'TMC-1'}, 7, 'AB', True, False, False, # will fail due to partial sequence
            marks=pytest.mark.xfail(
                raises=PartialBlockSequence,
                reason='Partial sequence repeat needed to get odd number block out of AB, but partial blocks are disabled',
                strict=True,
            )
        ),
        pytest.param(
            'monogrp_mpd_tmc', {'head':'MPD-1', 'tail' : 'TMC-1'}, 12, 'ABC', True, True, False, # will fail due to 3-functional TMC middle monomer as C
            marks=pytest.mark.xfail(
                raises=MorphologyError,
                reason='One of the monomers requested is non-linear (3-functional)',
                strict=True,
            )
        ),
        # PEG-PLGA
        ('monogrp_peg_plga', {}, 15, 'ABC', True, True, False), # test autogen
        ('monogrp_peg_plga', {}, 17, 'ABC', True, False, False), # test autogen with whole sequence
        ('monogrp_peg_plga', {'head':'PGA-1A', 'tail' : 'PGA-1B'}, 15, 'ABC', True, True, False), # test more complex sequence with non-default explicit end groups
        pytest.param(
            'monogrp_peg_plga', {'head':'PGA-1A', 'tail' : 'PGA-1B'}, 15, 'ABC', True, False, False, # will fail due to partial sequence
            marks=pytest.mark.xfail(
                raises=PartialBlockSequence,
                reason='Partial sequence repeat needed to get odd number block out of AB, but partial blocks are disabled',
                strict=True,
            )
        ),
        ('monogrp_peg_plga', {}, 40, 'ABCB', True, True, True), # test longer energy min
    ]
)
def test_build_linear_polymer_contributions(
    monomers : MonomerGroup,
    term_orient : dict[str, str],
    n_monomers : int,
    sequence : str,
    minimize_sequence : bool,
    allow_partial_sequences : bool,
    energy_minimize : bool,
    request : pytest.FixtureRequest, # allows for fixture expansion in parameterized arguments
) -> None:
    '''Test linear polymer builder assembles chains with right number of
    repeat units AND right number of atoms per repeat unit'''
    monomers = request.getfixturevalue(monomers) # unpack fixtures into their respective values
    monomers.term_orient = term_orient # this edit makes it VITAL that fixtures be function-level
    
    polymer = build_linear_polymer(
        monomers=monomers,
        n_monomers=n_monomers,
        sequence=sequence,
        minimize_sequence=minimize_sequence,
        allow_partial_sequences=allow_partial_sequences,
        add_Hs=False,
        energy_minimize=energy_minimize,
    )
    
    # characterize middle monomers
    n_rep_units = len(polymer.children)
    repeat_unit_sizes : dict[str, int] = {}
    repeat_unit_counts = Counter() # TODO: make use of this for checks!!
    for middle_monomers in polymer.children:
        repeat_unit_sizes[middle_monomers.name] = middle_monomers.n_particles
        repeat_unit_counts[middle_monomers.name] += 1
        
    # characterize end groups
    end_groups_requested = set(resname for head_or_tail, (resname, mol) in monomers.linear_end_groups().items())
    end_groups_used = set()
    for end_group in polymer.end_groups:
        if end_group is not None:
            end_groups_used.add(end_group.name)
            repeat_unit_sizes[end_group.name] = end_group.n_particles
            repeat_unit_counts[middle_monomers.name] += 1
    
    total_reps_match = (n_rep_units == n_monomers)
    contribs_match = all(num_monomers == monomers.contributions()[resname][0]
        for resname, num_monomers in repeat_unit_sizes.items()
    )
    end_groups_correct = (end_groups_used == end_groups_requested)
    # counts_match = ...
    
    assert all([total_reps_match, contribs_match, end_groups_correct]) #, and counts_match    )


@pytest.mark.parametrize(
    'sequence, term_orient, n_monomers, smiles_expected',
    compile_sequence_test_inputs(kernel_repeats=(2,))
)
def test_build_linear_polymer_sequence(
    monogrp_halogen_marked : MonomerGroup,
    sequence : str,
    term_orient : dict[str, str],
    n_monomers : int,
    smiles_expected : str,
) -> None:
    '''Test that repeat units are assembled in the expected order for a given input to the linear polymer builder'''
    mol_expected = explicit_mol_from_SMILES(smiles_expected)

    monogrp_halogen_marked.term_orient = term_orient
    chain = build_linear_polymer(
        monogrp_halogen_marked,
        n_monomers=n_monomers,
        sequence=sequence,
        allow_partial_sequences=False,
        add_Hs=False,
        energy_minimize=False, # No point minimzing, since we only care about order of groups, not coordinates
    )
    mol_actual = mbmol_to_rdmol(chain)

    assert mol_expected.HasSubstructMatch(mol_actual) and mol_actual.HasSubstructMatch(mol_expected) 
