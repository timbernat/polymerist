'''Unit tests for PDB file atom I/O utils'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import pytest

from typing import Any
from polymerist.molfiles.pdb.pdbatoms import SerialAtomLabeller, PDB_ATOM_RECORD_TOKENS, parse_pdb_atom_record


# TESTING ATOM READING
SAMPLE_PDB_RECORDS_WITH_INFO : dict[str, dict[str, Any]] = {
    # AA residue atoms
    'ATOM    189  C99 OCT     5      39.590  30.100  38.320  1.00  0.00' : {'Residue atom type': 'ATOM',
        'Atom serial number': 189,
        'Atom name': 'C99',
        'Alternate location indicator': '',
        'Residue name': 'OCT',
        'Chain identifier': '',
        'Residue sequence number': 5,
        'Residue insertion code': '',
        'X (angstrom)': 39.59,
        'Y (angstrom)': 30.1,
        'Z (angstrom)': 38.32,
        'Occupancy': 1.0,
        'Temperature factor': 0.0,
        'Segment identifier': '',
        'Element symbol': '',
        'Charge': '',
    },
    'ATOM    190 C100 OCT     5      38.850  31.110  37.700  1.00  0.00' : {
        'Residue atom type': 'ATOM',
        'Atom serial number': 190,
        'Atom name': 'C100',
        'Alternate location indicator': '',
        'Residue name': 'OCT',
        'Chain identifier': '',
        'Residue sequence number': 5,
        'Residue insertion code': '',
        'X (angstrom)': 38.85,
        'Y (angstrom)': 31.11,
        'Z (angstrom)': 37.7,
        'Occupancy': 1.0,
        'Temperature factor': 0.0,
        'Segment identifier': '',
        'Element symbol': '',
        'Charge': '',
    },
    'ATOM    246  OXT THR A  29      -0.317  20.109  12.824  1.00 25.00           O' : {
        'Residue atom type': 'ATOM',
        'Atom serial number': 246,
        'Atom name': 'OXT',
        'Alternate location indicator': '',
        'Residue name': 'THR',
        'Chain identifier': 'A',
        'Residue sequence number': 29,
        'Residue insertion code': '',
        'X (angstrom)': -0.317,
        'Y (angstrom)': 20.109,
        'Z (angstrom)': 12.824,
        'Occupancy': 1.0,
        'Temperature factor': 25.0,
        'Segment identifier': '',
        'Element symbol': 'O',
        'Charge': '',
    },
    'ATOM     17  N   GLN A   3      47.620  28.367   8.973  1.00 15.00           N' : {
        'Residue atom type': 'ATOM',
        'Atom serial number': 17,
        'Atom name': 'N',
        'Alternate location indicator': '',
        'Residue name': 'GLN',
        'Chain identifier': 'A',
        'Residue sequence number': 3,
        'Residue insertion code': '',
        'X (angstrom)': 47.62,
        'Y (angstrom)': 28.367,
        'Z (angstrom)': 8.973,
        'Occupancy': 1.0,
        'Temperature factor': 15.0,
        'Segment identifier': '',
        'Element symbol': 'N',
        'Charge': '',
    },
    ## testing altloc indicator
    'ATOM    890  OE1AGLN A1100      20.292   2.717  57.945  0.55 19.03           O' : { 
        'Residue atom type': 'ATOM',
        'Atom serial number': 890,
        'Atom name': 'OE1',
        'Alternate location indicator': 'A',
        'Residue name': 'GLN',
        'Chain identifier': 'A',
        'Residue sequence number': 1100,
        'Residue insertion code': '',
        'X (angstrom)': 20.292,
        'Y (angstrom)': 2.717,
        'Z (angstrom)': 57.945,
        'Occupancy': 0.55,
        'Temperature factor': 19.03,
        'Segment identifier': '',
        'Element symbol': 'O',
        'Charge': '',
    },
    ## testing residue insertion code
    'ATOM  11919  N   ARG D 100A     -9.676  74.726 -19.958  1.00105.71           N' : {  
        'Residue atom type': 'ATOM',
        'Atom serial number': 11919,
        'Atom name': 'N',
        'Alternate location indicator': '',
        'Residue name': 'ARG',
        'Chain identifier': 'D',
        'Residue sequence number': 100,
        'Residue insertion code': 'A',
        'X (angstrom)': -9.676,
        'Y (angstrom)': 74.726,
        'Z (angstrom)': -19.958,
        'Occupancy': 1.0,
        'Temperature factor': 105.71,
        'Segment identifier': '',
        'Element symbol': 'N',
        'Charge': '',
    },
    # non-AA residue atoms
    'HETATM   47  H21 UNL     1      49.668  24.248  10.436  1.00  0.00           H ' : {
        'Residue atom type': 'HETATM',
        'Atom serial number': 47,
        'Atom name': 'H21',
        'Alternate location indicator': '',
        'Residue name': 'UNL',
        'Chain identifier': '',
        'Residue sequence number': 1,
        'Residue insertion code': '',
        'X (angstrom)': 49.668,
        'Y (angstrom)': 24.248,
        'Z (angstrom)': 10.436,
        'Occupancy': 1.0,
        'Temperature factor': 0.0,
        'Segment identifier': '',
        'Element symbol': 'H',
        'Charge': '',
    },
    'HETATM 1071 FE   HEM A   1       8.128   7.371 -15.022 24.00 16.74          FE' : {
        'Residue atom type': 'HETATM',
        'Atom serial number': 1071,
        'Atom name': 'FE',
        'Alternate location indicator': '',
        'Residue name': 'HEM',
        'Chain identifier': 'A',
        'Residue sequence number': 1,
        'Residue insertion code': '',
        'X (angstrom)': 8.128,
        'Y (angstrom)': 7.371,
        'Z (angstrom)': -15.022,
        'Occupancy': 24.0,
        'Temperature factor': 16.74,
        'Segment identifier': '',
        'Element symbol': 'FE',
        'Charge': '',
    },
    ## testing charge sign
    'HETATM    1  C   LIG     1       0.000   0.000   0.000  1.00  0.00          C1+' : {
        'Residue atom type': 'HETATM',
        'Atom serial number': 1,
        'Atom name': 'C',
        'Alternate location indicator': '',
        'Residue name': 'LIG',
        'Chain identifier': '',
        'Residue sequence number': 1,
        'Residue insertion code': '',
        'X (angstrom)': 0.0,
        'Y (angstrom)': 0.0,
        'Z (angstrom)': 0.0,
        'Occupancy': 1.0,
        'Temperature factor': 0.0,
        'Segment identifier': '',
        'Element symbol': 'C1',
        'Charge': '+',
    }, 
    ## testing charge read when extraneous info beyond 80 lines is included
    'HETATM   15  C15 BAS A   1      10.098   3.266  -1.245  1.00  0.00           C -0.03007329' : {
        'Residue atom type': 'HETATM',
        'Atom serial number': 15,
        'Atom name': 'C15',
        'Alternate location indicator': '',
        'Residue name': 'BAS',
        'Chain identifier': 'A',
        'Residue sequence number': 1,
        'Residue insertion code': '',
        'X (angstrom)': 10.098,
        'Y (angstrom)': 3.266,
        'Z (angstrom)': -1.245,
        'Occupancy': 1.0,
        'Temperature factor': 0.0,
        'Segment identifier': '',
        'Element symbol': 'C',
        'Charge': '-',
    } 
}

@pytest.mark.parametrize('pdb_atom_record,expected_output', list(SAMPLE_PDB_RECORDS_WITH_INFO.items()))
def test_pdb_atom_record_parse_info_extraction(pdb_atom_record : str, expected_output : dict[str, Any]) -> None:
    '''Test that the pieces of information taken from PDB atom records take on the values expected'''
    assert parse_pdb_atom_record(pdb_atom_record) == expected_output

@pytest.mark.parametrize('pdb_atom_record', list(SAMPLE_PDB_RECORDS_WITH_INFO.keys()))
def test_pdb_atom_record_parse_type_coercion(pdb_atom_record : str) -> None:
    '''Test that the types assigned to pieces of information taken from PDB atom records are assigned the types expected'''
    assert all(
        type(field_value) == PDB_ATOM_RECORD_TOKENS[field_name][-1]
            for field_name, field_value in parse_pdb_atom_record(pdb_atom_record).items()
    )

# TESTING ATOM WRITING - LABELLING
ELEMS : tuple[str] = ('C', 'H', 'H', 'H', 'N', 'H', 'C', 'O', 'Cl') # atoms for methylcarbamoyl chloride (MCC)
@pytest.mark.parametrize(
    'mol_atom_elems, atom_label_width, include_elem_idx, default_elem_idx, expected_labels',
    [
        (ELEMS, 4, True , 0, ['C000', 'H000', 'H001', 'H002', 'N000', 'H003', 'C001', 'O000', 'Cl00']), # test with default PDD-compatible settings
        (ELEMS, 4, True , 2, ['C002', 'H002', 'H003', 'H004', 'N002', 'H005', 'C003', 'O002', 'Cl02']), # test element index offset
        (ELEMS, 3, True , 2, ['C02', 'H02', 'H03', 'H04', 'N02', 'H05', 'C03', 'O02', 'Cl2']), # test shorter atom label width
        (ELEMS, 1, True , 0, ['C', 'H', 'H', 'H', 'N', 'H', 'C', 'O', 'C']), # test truncation works below threshold where indices can be written
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
    '''Test that atom labelling behaves as expected with various label formatting configurations'''
    labeller = SerialAtomLabeller(
        atom_label_width=atom_label_width,
        include_elem_idx=include_elem_idx,
        default_elem_idx=default_elem_idx,
    )
    assert [labeller.get_atom_label(elem) for elem in mol_atom_elems] == expected_labels