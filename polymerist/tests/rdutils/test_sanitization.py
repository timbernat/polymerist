'''Unit tests for `sanitization` package'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import pytest

from typing import Generator
from rdkit import Chem

from polymerist.rdutils.sanitization import sanitize_mol_outputs


# helper functions
@sanitize_mol_outputs
def makemols_generator(smiles : list[str]) -> Generator[Chem.Mol, None, None]:
    for smi in smiles:
        yield Chem.MolFromSmiles(smi, sanitize=False)

@sanitize_mol_outputs
def makemols_container(smiles : list[str]) -> list[Chem.Mol]:
    ms = []
    for smi in smiles:
        ms.append(Chem.MolFromSmiles(smi, sanitize=False))
    return ms

@sanitize_mol_outputs
def makemols_singular(smiles : list[str]) -> list[Chem.Mol]:
    return Chem.MolFromSmiles(smiles[0])