'''Tests construction of structures for linear copolymers (and relevant subfamilies, e.g. homopolymers)'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import pytest

from collections import Counter

from polymerist.polymers.monomers.repr import MonomerGroup
from polymerist.polymers.monomers.fragments import PE_FRAGMENTS, MPD_TMC_FRAGMENTS, PEG_PLGA_FRAGMENTS

from polymerist.polymers.building import build_linear_polymer
from polymerist.polymers.exceptions import MorphologyError, PartialBlockSequence, EmptyBlockSequence


@pytest.fixture(scope='function')
def monogrp_polyethylene() ->  MonomerGroup:
    return MonomerGroup(monomers=PE_FRAGMENTS)

@pytest.fixture(scope='function')
def monogrp_mpd_tmc() ->  MonomerGroup:
    return MonomerGroup(monomers=MPD_TMC_FRAGMENTS)

@pytest.fixture(scope='function')
def monogrp_peg_plga() ->  MonomerGroup:
    return MonomerGroup(monomers=PEG_PLGA_FRAGMENTS)


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
def test_build_linear_polymer(
        monomers : MonomerGroup,
        term_orient : dict[str, str],
        n_monomers : int,
        sequence : str,
        minimize_sequence : bool,
        allow_partial_sequences : bool,
        energy_minimize : bool,
        request : pytest.FixtureRequest, # allows for fixture expansion in parameterized arguments
    ) -> None:
    '''Test linear polymer builder behavior under varing sets of parameters'''
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
    residue_sizes : dict[str, int] = {}
    residue_counts = Counter() # TODO: make use of this for checks!!
    for middle_monomers in polymer.children:
        residue_sizes[middle_monomers.name] = middle_monomers.n_particles
        residue_counts[middle_monomers.name] += 1
        
    # characterize end groups
    end_groups_requested = set(resname for head_or_tail, (resname, mol) in monomers.linear_end_groups().items())
    end_groups_used = set()
    for end_group in polymer.end_groups:
        if end_group is not None:
            end_groups_used.add(end_group.name)
            residue_sizes[end_group.name] = end_group.n_particles
            residue_counts[middle_monomers.name] += 1
    
    total_reps_match = (n_rep_units == n_monomers)
    contribs_match = all(num_monomers == monomers.contributions()[resname][0]
        for resname, num_monomers in residue_sizes.items()
    )
    end_groups_correct = (end_groups_used == end_groups_requested)
    # counts_match = ...
    
    assert all([total_reps_match, contribs_match, end_groups_correct]) #, and counts_match    )