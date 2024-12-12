'''Unit tests for PDB file I/O utils'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import pytest
from polymerist.molfiles.pdb import SerialAtomLabeller


ELEMS : tuple[str] = ('C', 'H', 'H', 'H', 'N', 'H', 'C', 'O', 'Cl') # atoms for methylcarbamoyl chloride (MCC)

@pytest.mark.parametrize(
    'mol_atom_elems, atom_label_width, include_elem_idx, default_elem_idx, expected_labels',
    [
        (ELEMS, 4, True, 0, ['C000', 'H000', 'H001', 'H002', 'N000', 'H003', 'C001', 'O000', 'Cl00']), # test with default PDD-compatible settings
        (ELEMS, 4, True, 2, ['C002', 'H002', 'H003', 'H004', 'N002', 'H005', 'C003', 'O002', 'Cl02']), # test element index offset
        (ELEMS, 3, True, 2, ['C02', 'H02', 'H03', 'H04', 'N02', 'H05', 'C03', 'O02', 'Cl2']), # test shorter atom label width
        (ELEMS, 1, True, 0, ['C', 'H', 'H', 'H', 'N', 'H', 'C', 'O', 'C']), # test truncation works below threshold where indices can be written
        (ELEMS, 4, False, 0, ['C   ', 'H   ', 'H   ', 'H   ', 'N   ', 'H   ', 'C   ', 'O   ', 'Cl  ']), # test without element indices
        (ELEMS, 4, False, 7, ['C   ', 'H   ', 'H   ', 'H   ', 'N   ', 'H   ', 'C   ', 'O   ', 'Cl  ']), # test that default indices has no impact when indices aren't present
        (ELEMS, 0, False, 0, ['', '', '', '', '', '', '', '', '']), # test null-width labels
        # Invalid input handling checks
        pytest.param(
            ELEMS, -1, True, 0, [], # test that negative label width is rejected as intended
            marks=pytest.mark.xfail(
                raises=ValueError,
                reason='Negative atom label widths not allowed',
                strict=True,
            )
        ),
        pytest.param(
            ELEMS, 4, True, -5, [], # test that negative default indices are rejected as intended
            marks=pytest.mark.xfail(
                raises=ValueError,
                reason='Negative element indices not allowed',
                strict=True,
            )
        ),
        pytest.param(
            tuple(len(elem) for elem in ELEMS), 4, True, 0, [], # test that negative default indices are rejected as intended
            marks=pytest.mark.xfail(
                raises=TypeError,
                reason='Must pass atom elements as strings',
                strict=True,
            )
        ),
    ]
)
def test_atom_labeller(
        mol_atom_elems : tuple[str],
        atom_label_width : int,
        include_elem_idx : bool,
        default_elem_idx : int,
        expected_labels : list[str],
    ) -> None:
    '''Test that atom labelling hebaves as expected with various label formatting configurations'''
    labeller = SerialAtomLabeller(
        atom_label_width=atom_label_width,
        include_elem_idx=include_elem_idx,
        default_elem_idx=default_elem_idx,
    )
    assert [labeller.get_atom_label(elem) for elem in mol_atom_elems] == expected_labels