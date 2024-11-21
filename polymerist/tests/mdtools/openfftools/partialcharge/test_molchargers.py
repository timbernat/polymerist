'''Unit tests for `molchargers` package'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import pytest

from openff.toolkit import Molecule

from polymerist.polymers.monomers.specification import expanded_SMILES
from polymerist.mdtools.openfftools.partialcharge.molchargers import MolCharger
from polymerist.mdtools.openfftools.partialcharge.rescharge.interface import LibraryCharger


# Test MolCharger subclass registration
def test_molcharger_registers_subclasses() -> None:
    '''Test that the MolCharger tracks subclasses'''
    assert hasattr(MolCharger, 'subclass_registry')
    
@pytest.mark.parametrize('expected_charge_method_name, molcharger_subclass', MolCharger.subclass_registry.items()) # NOTE: this will fail if test_molcharger_registers_subclasses() fails
def test_molcharger_subclass_attr_registration(molcharger_subclass : type[MolCharger], expected_charge_method_name : str) -> None:
    '''Test that all MolCharger subclasses define and are registered under their "CHARGING_METHOD" class property'''
    assert hasattr(molcharger_subclass, 'CHARGING_METHOD') and (getattr(molcharger_subclass, 'CHARGING_METHOD') == expected_charge_method_name)


# Test MolCharger subclass implementations
@pytest.fixture
def offmol() -> Molecule:
    '''Dummy Molecule object for testing'''
    # DEV: worthing double-checking that partial charges are initially empty??
    return Molecule.from_smiles('c1ccccc1C(=O)O') # benzoic acid - nice and small, but with some non-trivial structure
    
MOLCHARGER_TYPES_TO_TEST = [
    molcharger_subclass
        for molcharger_subclass in MolCharger.subclass_registry.values() # LibraryCharger behave differently to other MolCharger and are kind of a pain in the ass generally...
            if molcharger_subclass != LibraryCharger                     # ...intend to deprecate and revamp them eventually, so will just exclude them from testing for now
]
    
@pytest.mark.parametrize('molcharger_subclass', MOLCHARGER_TYPES_TO_TEST)
def test_molchargers_assign_charges(offmol : Molecule, molcharger_subclass : type[MolCharger]) -> None:
    charger = molcharger_subclass()
    cmol = charger.charge_molecule(offmol)
    assert cmol.partial_charges is not None # should assign charges to the new, copied molecule

@pytest.mark.parametrize('molcharger_subclass', MOLCHARGER_TYPES_TO_TEST)
def test_molchargers_act_readonly(offmol : Molecule, molcharger_subclass : type[MolCharger]) -> None:
    charger = molcharger_subclass()
    cmol = charger.charge_molecule(offmol)
    assert offmol.partial_charges is None # should NOT affect the 

@pytest.mark.parametrize('molcharger_subclass', MOLCHARGER_TYPES_TO_TEST)
def test_molchargers_record_charge_method(offmol : Molecule, molcharger_subclass : type[MolCharger]) -> None:
    charger = molcharger_subclass()
    cmol = charger.charge_molecule(offmol)
    
    recorded_charge_method = cmol.properties.get('charge_method', None)
    assert (recorded_charge_method is not None) and (recorded_charge_method == getattr(molcharger_subclass, 'CHARGING_METHOD'))
