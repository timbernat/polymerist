'''Unit tests for `molchargers` package'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import pytest

# NOTE: inverted custom imports here to get polymerist.openfftools import first, 
# Done so that if OpenFF is not found, a helpful installation error with be raised prior to attempting direct openff.toolkit imports below
from polymerist.genutils.importutils.dependencies import modules_installed
from polymerist.mdtools.openfftools.partialcharge.molchargers import MolCharger, ABE10Charger, EspalomaCharger, NAGLCharger
from polymerist.mdtools.openfftools.partialcharge.rescharge.interface import LibraryCharger

from openff.toolkit import Molecule
from openff.toolkit.utils.toolkits import OPENEYE_AVAILABLE


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
    
## selectively register test to avoid failures due to missing optional dependencies
MOLCHARGER_TYPES_TO_TEST : list[type[MolCharger]] = [] 
if modules_installed('openeye.oechem', 'openeye.oeomega') and OPENEYE_AVAILABLE: # extra check needed block check when missing license (as is the case for the open-source polymerist repo)
    MOLCHARGER_TYPES_TO_TEST.append(ABE10Charger)
if modules_installed('espaloma_charge'):
    MOLCHARGER_TYPES_TO_TEST.append(EspalomaCharger)
if modules_installed('openff.nagl', 'openff.nagl_models'):
    MOLCHARGER_TYPES_TO_TEST.append(NAGLCharger)
# MOLCHARGER_TYPES_TO_TEST.append(
    # LibraryCharger  # LibraryCharger behave differently to other MolCharger and are kind of a pain in the ass generally...
# )                   # ...intend to deprecate and revamp them eventually, so will just exclude them from testing for now
    
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
