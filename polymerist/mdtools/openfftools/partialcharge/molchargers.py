'''Classes for partial charge assignment of OpenFF Molecules'''

import logging
LOGGER = logging.getLogger(__name__)

from typing import Any, ClassVar, Union
from abc import ABC, abstractmethod, abstractproperty

from rdkit import Chem
from openff.units import unit as offunit
from openff.toolkit.topology.molecule import Molecule
from openff.toolkit.utils.exceptions import ToolkitUnavailableException # TODO : use chargemethods.TOOLKITS_BY_CHARGE_METHOD to automatically determine whether/which toolkits are available for each method

from .. import TKREGS, _OE_TKWRAPPER_IS_AVAILABLE, OEUnavailableException
from .chargemethods import NAGL_MODEL
from ....genutils.decorators.classmod import register_subclasses
from ....genutils.decorators.functional import optional_in_place


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
        if not _OE_TKWRAPPER_IS_AVAILABLE:
            raise OEUnavailableException # AM1-BCC-ELF10 is exclusively available thru OpenEye; if it is not present, then, must err
        uncharged_mol.assign_partial_charges(partial_charge_method='am1bccelf10', toolkit_registry=TKREGS['OpenEye Toolkit']) # TODO : provide support for AMBER / RDKit if OE license is unavailable

class EspalomaCharger(MolCharger):
    '''Charger class for EspalomaCharge charging'''
    CHARGING_METHOD : ClassVar[str] = 'Espaloma-AM1-BCC'

    @optional_in_place
    def _charge_molecule(self, uncharged_mol : Molecule) -> None:
        uncharged_mol.assign_partial_charges(partial_charge_method='espaloma-am1bcc', toolkit_registry=TKREGS['Espaloma Charge Toolkit'])

class NAGLCharger(MolCharger):
    '''Charger class for NAGL charging'''
    CHARGING_METHOD : ClassVar[str] = 'NAGL'

    @optional_in_place
    def _charge_molecule(self, uncharged_mol : Molecule) -> None:
        nagl_charges = NAGL_MODEL.compute_property(uncharged_mol, check_domains=True, error_if_unsupported=True)
        uncharged_mol.partial_charges = nagl_charges * offunit.elementary_charge # need to have OpenFF-style units attached to set "partial_charges" property
