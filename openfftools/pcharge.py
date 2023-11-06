'''Custom typehints and classes for partial charge assignment of OpenFF Molecules'''

import logging
LOGGER = logging.getLogger(__name__)

from typing import Any, ClassVar
from abc import ABC, abstractmethod, abstractproperty
from openff.toolkit.topology.molecule import Molecule

from . import TKREGS
from ..genutils.decorators.classmod import register_subclasses
from ..genutils.decorators.functional import optional_in_place


# ABSTRACT AND CONCRETE CLASSES FOR CHARGING MOLECULES
@register_subclasses(key_attr='CHARGING_METHOD')
class MolCharger(ABC):
    '''Base interface for defining various methods of generating and storing atomic partial charges'''
    @abstractproperty
    @classmethod
    def CHARGING_METHOD(cls):
        '''For setting the name of the method as a class attribute in child classes'''
        pass

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
class ABE10Charger(MolCharger):
    '''Charger class for AM1-BCC-ELF10 exact charging'''
    CHARGING_METHOD : ClassVar[str] = 'AM1-BCC-ELF10'

    @optional_in_place
    def _charge_molecule(self, uncharged_mol : Molecule) -> None:
        '''Concrete implementation for AM1-BCC-ELF10'''
        uncharged_mol.assign_partial_charges(partial_charge_method='am1bccelf10', toolkit_registry=TKREGS['OpenEye Toolkit'])

class EspalomaCharger(MolCharger):
    '''Charger class for EspalomaCharge charging'''
    CHARGING_METHOD : ClassVar[str] = 'Espaloma-AM1-BCC'

    @optional_in_place
    def _charge_molecule(self, uncharged_mol : Molecule) -> None:
        uncharged_mol.assign_partial_charges(partial_charge_method='espaloma-am1bcc', toolkit_registry=TKREGS['Espaloma Charge Toolkit'])