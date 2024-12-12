'''Testing that copolymer sequencing scales (and fails) as expected'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Any
from dataclasses import asdict

import pytest
from pathlib import Path

from polymerist.polymers.building.sequencing import LinearCopolymerSequencer as LCS
from polymerist.polymers.exceptions import EmptyBlockSequence, PartialBlockSequence, InsufficientChainLength, EndGroupDominatedChain


@pytest.fixture
def sequencer() -> LCS:
    '''A sample sequencer with known, valid inputs'''
    return LCS(sequence_kernel='ABAB', n_repeat_units=14, n_repeat_units_terminal=2)

@pytest.mark.parametrize(
    'inputs',
    [
        {
            'sequence_kernel' : 'AB',
            'n_repeat_units' : 10,
            'n_repeat_units_terminal' : 1
        },
        pytest.param(
            {
                'sequence_kernel' : 'BAC',
                'n_repeat_units' : 1,
                'n_repeat_units_terminal' : 2
            },
            marks=pytest.mark.xfail(
                raises=EndGroupDominatedChain,
                reason='Results in (unsatisfiable) negative number of middle monomers',
                strict=True,
            )
        ),
        pytest.param(
            {
                'sequence_kernel' : '',
                'n_repeat_units' : 7,
                'n_repeat_units_terminal' : 1
            },
            marks=pytest.mark.xfail(
                raises=EmptyBlockSequence,
                reason='No sequence kernel provided',
                strict=True,
            )
        ), 
    ]
)
def test_LCS_input_validation(inputs : dict[str, Any]) -> None:
    '''Test that invalid Sequence input are correctly rejected'''
    _ = LCS(**inputs) # no assert needed, just checking when initialization completes
    
def test_LCS_copying(sequencer : LCS) -> None:
    '''Test that sequencers are properly copied in a read-only manner'''
    sequencer_clone = sequencer.copy()
    
    # tamper with the parameters of the copy in a way that guarantees distinctness
    sequencer_clone.sequence_kernel = 2*sequencer.sequence_kernel
    sequencer_clone.n_repeat_units += 2
    sequencer_clone.n_repeat_units_terminal += 1
    
    # check that the original WASN'T tampered with
    assert asdict(sequencer) != asdict(sequencer_clone)
    
    
@pytest.mark.parametrize(
    'sequencer, expected_kernel',
    [
        (LCS('ABC', n_repeat_units=12), 'ABC') , # test irrreducible case
        (LCS('ABAB', n_repeat_units=12), 'AB'), # test unreduced case
    ]
)
def test_LCS_reduction(sequencer : LCS, expected_kernel : str) -> None:
    '''Test that shortest repeating subsequences of sequencer kernels are correctly identified'''
    sequencer.reduce()
    assert sequencer.sequence_kernel == expected_kernel
    
@pytest.mark.parametrize(
    'sequencer, allow_partials, expected_sequence, expected_length',
    [
        # tests for homopolymers
        (LCS('A', 5, 1), True , 'A', 4),
        (LCS('A', 5, 1), False, 'A', 4), # partial block single-monomer sequence will never exist, so "allow_partial_sequences" setting shouldn't matter)
        pytest.param(
            LCS('A', 1, 1), True, 'A', 1, # test that all-end group (i.e. no middle monomer) case is correctly rejected
            marks=pytest.mark.xfail(
                raises=InsufficientChainLength,
                reason='No middle monomers can be accomodated',
                strict=True,
            ),
        ),
        # tests for "true" copolymers
        (LCS('ABC', 10, 2), True, 'ABCABCAB', 1),
        pytest.param(
            LCS('ABC', 10, 2), False, 'ABCABCAB', 1, # test that partial-sequence ban correctly blocks partial sequences...
            marks=pytest.mark.xfail(
                raises=PartialBlockSequence,
                reason='Partial sequence repeats have not been allowed',
                strict=True,
            ),
        ),
        (LCS('ABC', 11, 2), False, 'ABC', 3), # ...unless the resulting sequence happens to be a whole multiple
        pytest.param(
            LCS('ABC', 2, 2), True, '', 1, # test that all-end group (i.e. no middle monomer) case is correctly rejected...
            marks=pytest.mark.xfail(
                raises=InsufficientChainLength,
                reason='No middle monomers can be accomodated',
                strict=True,
            ),
        ),
        (LCS('ABC', 4, 2), True, 'AB', 1), # ... and finally, check that nonempty sequences SMALLER than the kernel are also recognized if partials are permitted
    ]
)
def test_LCS_procrustean_alignment(sequencer : LCS, allow_partials : bool, expected_sequence : str, expected_length : int) -> None:
    '''Test capability (and prechecks) for fitting sequence to target chain length'''
    seq, n_reps = sequencer.procrustean_alignment(allow_partial_sequences=allow_partials)
    assert (seq == expected_sequence) and (n_reps == expected_length)