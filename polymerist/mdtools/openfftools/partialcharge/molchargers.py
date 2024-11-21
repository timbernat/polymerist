'''Classes for partial charge assignment of OpenFF Molecules'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import logging
LOGGER = logging.getLogger(__name__)

from typing import Union
from abc import ABC, abstractmethod

from rdkit import Chem
from openff.toolkit.topology.molecule import Molecule
from openff.toolkit.utils.exceptions import ToolkitUnavailableException

from ....genutils.importutils.dependencies import requires_modules
from ....genutils.decorators.functional import optional_in_place
from ....genutils.decorators.classmod import register_subclasses, register_abstract_class_attrs


def has_partial_charges(mol : Union[Molecule, Chem.Mol]) -> bool:
    '''Check if a molecular representation (either a OpenFF Molecule or and RDKit Mol) has partial charges assigned'''
    if isinstance(mol, Molecule):
        return (mol.partial_charges is not None)
    elif isinstance(mol, Chem.Mol): # NOTE : would work just as well with just an "if" here, but elif communicates intent better
        return bool(mol.HasProp('atom.dprop.PartialCharge')) # RDKit returns as int 0 or 1 instead of bool for some reason
    else:
        raise TypeError(f'Cannot check partial charge status of object of type {mol.__class__.__name__}')
    
    
# ABSTRACT AND CONCRETE CLASSES FOR CHARGING MOLECULES
@register_subclasses(key_attr='CHARGING_METHOD')
@register_abstract_class_attrs('CHARGING_METHOD')
class MolCharger(ABC):
    '''Base interface for defining various methods of generating and storing atomic partial charges'''
    @abstractmethod
    @optional_in_place
    def _charge_molecule(self, uncharged_mol : Molecule) -> None: 
        '''Method for assigning molecular partial charges - concrete implementation in child classes'''
        pass

    @optional_in_place
    def charge_molecule(self, uncharged_mol : Molecule) -> None:
        '''Wraps charge method call with logging'''
        LOGGER.info(f'Assigning partial charges via the "{self.CHARGING_METHOD}" method')
        self._charge_molecule(uncharged_mol, in_place=True) # must be called in-place for external optional_in_place wrapper to work as expected
        uncharged_mol.properties['charge_method'] = self.CHARGING_METHOD # record method within Molecule for reference
        LOGGER.info(f'Successfully assigned "{self.CHARGING_METHOD}" charges')

# CONCRETE IMPLEMENTATIONS OF DIFFERENT CHARGING METHODS
class ABE10Charger(MolCharger, CHARGING_METHOD= 'AM1-BCC-ELF10'):
    '''Charger class for AM1-BCC-ELF10 exact charging'''
    @requires_modules('openeye.oechem', 'openeye.oeomega', missing_module_error=ToolkitUnavailableException) # for whatever weird reason the toplevel openeye package has no module spec, so just checking "openeye" isn't enough
    @optional_in_place
    def _charge_molecule(self, uncharged_mol : Molecule) -> None:
        from openff.toolkit.utils.openeye_wrapper import OpenEyeToolkitWrapper
        
        uncharged_mol.assign_partial_charges(
            partial_charge_method='am1bccelf10',
            toolkit_registry=OpenEyeToolkitWrapper(), # instance init will raise exception if license or OpenEye packages are missing
        ) # TODO : find decent alternative if OpenEye license is missing (AmberTools doesn't do ELF10 and doesn't work on Windows)

class EspalomaCharger(MolCharger, CHARGING_METHOD='Espaloma-AM1-BCC'):
    '''Charger class for EspalomaCharge charging'''
    @requires_modules('espaloma_charge', missing_module_error=ToolkitUnavailableException)
    @optional_in_place
    def _charge_molecule(self, uncharged_mol : Molecule) -> None:
        from espaloma_charge.openff_wrapper import EspalomaChargeToolkitWrapper
        
        uncharged_mol.assign_partial_charges(
            partial_charge_method='espaloma-am1bcc', # this is actually the ONLY charge method the EspalomaChargeToolkitWrapper supports
            toolkit_registry=EspalomaChargeToolkitWrapper(),
        )

class NAGLCharger(MolCharger, CHARGING_METHOD='NAGL'):
    '''Charger class for NAGL charging'''
    @requires_modules('openff.nagl', 'openff.nagl_models', missing_module_error=ToolkitUnavailableException)
    @optional_in_place
    def _charge_molecule(self, uncharged_mol : Molecule) -> None:
        from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
        
        uncharged_mol.assign_partial_charges(
            partial_charge_method='openff-gnn-am1bcc-0.1.0-rc.3.pt', # 'openff-gnn-am1bcc-0.1.0-rc.2.pt',
            toolkit_registry=NAGLToolkitWrapper(),
        )
