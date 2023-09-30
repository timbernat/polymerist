'''Simplifies creation of Simulations which correspond to a particular thermodynamic ensembles'''

import logging
LOGGER = logging.getLogger(__name__)

from typing import Any, ClassVar, Iterable, Optional
from abc import ABC, abstractmethod, abstractproperty
from dataclasses import dataclass, field

from openmm import Integrator, VerletIntegrator, LangevinMiddleIntegrator
from openmm.openmm import Force, MonteCarloBarostat
from openmm.unit import Quantity, kelvin, atmosphere, picosecond

from ..genutils.fileutils import jsonio
from ..genutils.decorators.classmod import register_subclasses


# PARAMETER CLASSES
@dataclass
class ThermoParameters(jsonio.JSONifiable):
    '''For recording temperature, pressure, ensemble, and other thermodynamic parameters'''
    ensemble       : str = 'NVT'
    temperature    : Quantity = field(default_factory=lambda : (300 * kelvin)) # just specifying Quantities as default doesn't cut it, since these are (evidently) mutable defaults which dataclasses can't digest
    pressure       : Quantity = field(default_factory=lambda : (1 * atmosphere))

    friction_coeff : Quantity = field(default_factory=lambda : (1 / picosecond))
    barostat_freq  : int = 25

    def __post_init__(self) -> None:
        self.ensemble = self.ensemble.upper() # ensure ensemble name is upper-case

    # JSON serialization
    @staticmethod
    def serialize_json_dict(unser_jdict : dict[Any, Any]) -> dict[str, jsonio.JSONSerializable]:
        '''Serialize unit-ful Quantity attrs in a way the JSON can digest'''
        return jsonio.json_serialize_quantities(unser_jdict)
    
    @staticmethod
    def unserialize_json_dict(ser_jdict : dict[str, jsonio.JSONSerializable]) -> dict[Any, Any]:
        '''For unserializing unit-ful Quantities upon load from json file'''
        return jsonio.json_unserialize_quantities(ser_jdict)
    

# ABSTRACT BASE FOR CREATING ENSEMBLE-SPECIFIC SIMULATION
@dataclass
@register_subclasses(key_attr='ensemble')
class EnsembleFactory(ABC):
    '''Base class for implementing interface for generating ensemble-specific simulations'''
    thermo_params : ThermoParameters

    @classmethod
    def from_thermo_params(cls, thermo_params : ThermoParameters) -> 'EnsembleFactory':
        '''class method to automatically perform registry lookup (simplifies imports and use cases)'''
        return EnsembleFactory.subclass_registry[thermo_params.ensemble](thermo_params)

    # ABSTRACT METHODS AND PROPERTIES (TO BE IMPLEMENTED FOR EACH PARTICULAR ENSEMBLE)
    @abstractproperty
    @classmethod
    def ensemble(self) -> str:
        '''Specify state variables of ensemble'''
        pass

    @abstractproperty
    @classmethod
    def ensemble_name(self) -> str:
        '''Specify name of ensemble'''
        pass

    @abstractmethod
    def _integrator(self, time_step : Quantity) -> Integrator:
        '''Specify how to integrate forces in each timestep'''
        pass

    @abstractmethod
    def _forces(self) -> Optional[Iterable[Force]]:
        '''Specify any additional force contributions to position/velocity updates'''
        pass

    # CONCRETE METHODS WITH INTEGRATED LOGGING
    def integrator(self, time_step : Quantity) -> Integrator:
        '''Specify how to integrate forces in each timestep'''
        integrator = self._integrator(time_step)
        LOGGER.info(f'Created {integrator.__class__.__name__} for {self.desc}')
        
        return integrator

    def forces(self) -> Optional[Iterable[Force]]:
        '''Specify any additional force contributions to position/velocity updates'''
        forces = self._forces()
        if forces:
            force_str = ', '.join(force.__class__.__name__ for force in forces)
            LOGGER.info(f'Created {force_str} Force(s) for {self.desc}')

        return forces

    # REPRESENTATION AND PRETTY-PRINTING
    _REPR_ATTRS = ('ensemble', 'ensemble_name')

    @property
    def desc(self) -> str:
        '''Verbal description of ensemble'''
        return f'{self.ensemble} ({self.ensemble_name.capitalize()}) ensemble'

    def __repr__(self) -> str:
        '''Provide a description of the ensemble and mechanics used'''
        
        attr_str = ', '.join(
            f'{attr_name}={getattr(self, attr_name)}'    
                for attr_name in self._REPR_ATTRS
        )
        return f'{self.__class__.__name__}({attr_str})'
    

# CONCRETE IMPLEMENTATIONS OF ENSEMBLES
@dataclass
class NVESimulationFactory(EnsembleFactory):
    thermo_params : ThermoParameters

    ensemble : ClassVar[str] = 'NVE'
    ensemble_name : ClassVar[str] = 'microcanonical'

    def _integrator(self, time_step : Quantity) -> Integrator:
        return VerletIntegrator(stepSize=time_step)
    
    def _forces(self) -> Optional[Iterable[Force]]:
        return None
    
@dataclass
class NVTSimulationFactory(EnsembleFactory): # TODO : add implementation support for Andersen and Nose-Hoover thermostats (added to forces instead)
    thermo_params : ThermoParameters

    ensemble : ClassVar[str] = 'NVT'
    ensemble_name : ClassVar[str] = 'canonical'

    def _integrator(self, time_step : Quantity) -> Integrator:
        return LangevinMiddleIntegrator(self.thermo_params.temperature, self.thermo_params.friction_coeff, time_step)
    
    def _forces(self) -> Optional[Iterable[Force]]:
        return None
    
@dataclass
class NPTSimulationFactory(EnsembleFactory):
    thermo_params : ThermoParameters

    ensemble : ClassVar[str] = 'NPT'
    ensemble_name : ClassVar[str] = 'isothermal-isobaric'

    def _integrator(self, time_step : Quantity) -> Integrator:
        return LangevinMiddleIntegrator(self.thermo_params.temperature, self.thermo_params.friction_coeff, time_step)
    
    def _forces(self) -> Optional[Iterable[Force]]:
        return [MonteCarloBarostat(self.thermo_params.pressure, self.thermo_params.temperature, self.thermo_params.barostat_freq)]
    
  