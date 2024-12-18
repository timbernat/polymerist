'''Simplifies creation of Simulations which correspond to a particular thermodynamic ensembles'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import logging
LOGGER = logging.getLogger(__name__)

from typing import Any, ClassVar, Iterable, Optional
from abc import ABC, abstractmethod, abstractproperty
from dataclasses import dataclass, field

from openmm import Integrator, VerletIntegrator, LangevinMiddleIntegrator
from openmm.openmm import Force, MonteCarloBarostat
from openmm.unit import Quantity, kelvin, atmosphere, picosecond

from ...genutils.decorators.classmod import register_subclasses, register_abstract_class_attrs
from ...genutils.fileutils.jsonio.jsonify import make_jsonifiable
from ...genutils.fileutils.jsonio.serialize import QuantitySerializer


# PARAMETER CLASSES
@make_jsonifiable(type_serializer=QuantitySerializer)
@dataclass
class ThermoParameters:
    '''For recording temperature, pressure, ensemble, and other thermodynamic parameters'''
    ensemble       : str = 'NVT'
    temperature    : Quantity = field(default_factory=lambda : (300 * kelvin)) # just specifying Quantities as default doesn't cut it, since these are (evidently) mutable defaults which dataclasses can't digest
    pressure       : Quantity = field(default_factory=lambda : (1 * atmosphere))

    friction_coeff : Quantity = field(default_factory=lambda : (1 / picosecond))
    barostat_freq  : int = 25

    def __post_init__(self) -> None:
        self.ensemble = self.ensemble.upper() # ensure ensemble name is upper-case


# ABSTRACT BASE FOR CREATING ENSEMBLE-SPECIFIC SIMULATION
@dataclass
@register_subclasses(key_attr='ensemble')
@register_abstract_class_attrs('ensemble', 'ensemble_name')
class EnsembleFactory(ABC):
    '''Base class for implementing interface for generating ensemble-specific simulations'''
    thermo_params : ThermoParameters

    @classmethod
    def from_thermo_params(cls, thermo_params : ThermoParameters) -> 'EnsembleFactory':
        '''class method to automatically perform registry lookup (simplifies imports and use cases)'''
        return EnsembleFactory.subclass_registry[thermo_params.ensemble](thermo_params)

    # ENSEMBLE NAMING ATTRIBUTES
    @property
    def desc(self) -> str:
        '''Verbal description of ensemble'''
        return f'{self.ensemble} ({self.ensemble_name.capitalize()}) ensemble'

    def __post_init__(self) -> None:
        self.__doc__ = f'Factory class for the {self.desc}' # auto-generate docstring by ensemble

    # DEFINING ENSEMBLE-SPECIFIC INTEGRATORS AND FORCES
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


# CONCRETE IMPLEMENTATIONS OF ENSEMBLES
@dataclass
class NVESimulationFactory(EnsembleFactory, ensemble='NVE', ensemble_name='microcanonical'):
    def _integrator(self, time_step : Quantity) -> Integrator:
        return VerletIntegrator(time_step)
    
    def _forces(self) -> Optional[Iterable[Force]]:
        return None
    
@dataclass
class NVTSimulationFactory(EnsembleFactory, ensemble='NVT', ensemble_name='canonical'): 
    def _integrator(self, time_step : Quantity) -> Integrator:
        return LangevinMiddleIntegrator(self.thermo_params.temperature, self.thermo_params.friction_coeff, time_step)
    
    def _forces(self) -> Optional[Iterable[Force]]:
        return None # TODO : add implementation support for Andersen and Nose-Hoover thermostats (added to forces instead)
    
@dataclass
class NPTSimulationFactory(EnsembleFactory, ensemble='NPT', ensemble_name='isothermal-isobaric'):
    def _integrator(self, time_step : Quantity) -> Integrator:
        return LangevinMiddleIntegrator(self.thermo_params.temperature, self.thermo_params.friction_coeff, time_step)
    
    def _forces(self) -> Optional[Iterable[Force]]:
        return [MonteCarloBarostat(self.thermo_params.pressure, self.thermo_params.temperature, self.thermo_params.barostat_freq)]
    
  