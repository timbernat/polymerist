'''Simplifies creation of Simulations which correspond to a particular thermodynamic ensembles'''

import logging
LOGGER = logging.getLogger(__name__)

from abc import ABC, abstractmethod, abstractproperty
from typing import Any, Iterable, Optional, Union

from openff.units import unit as offunit # need OpenFF version of unit for Interchange positions for some reason
from openff.interchange import Interchange

from openmm import Integrator, VerletIntegrator, LangevinMiddleIntegrator
from openmm.openmm import Force, MonteCarloBarostat
from openmm.app import Simulation

from .records import SimulationParameters
from ..genutils.decorators.classmod import register_subclasses


# ABSTRACT BASE FOR CREATING ENSEMBLE-SPECIFIC SIMULATION
@register_subclasses(key_attr='ensemble')
class EnsembleSimulationFactory(ABC):
    '''Base class for implementing interface for generating ensemble-specific simulations'''
    # Abstract methods and properties
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
    def integrator(self, sim_params : SimulationParameters) -> Integrator:
        '''Specify how to integrate forces in each timestep'''
        pass

    @abstractmethod
    def forces(self, sim_params : SimulationParameters) -> Optional[Iterable[Force]]:
        '''Specify any additional force contributions to position/velocity updates'''
        pass

    # Concrete methods and properties
    _REPR_ATTRS = ('ensemble', 'ensemble_name')
    def __repr__(self) -> str:
        '''Provide a description of the ensemble and mechanics used'''
        
        attr_str = ', '.join(
            f'{attr_name}={getattr(self, attr_name)}'    
                for attr_name in self._REPR_ATTRS
        )
        return f'{self.__class__.__name__}({attr_str})'
    
    @property
    def desc(self) -> str:
        '''Verbal description of ensemble'''
        return f'{self.ensemble} ({self.ensemble_name.capitalize()} ensemble)'

    def create_simulation(self, interchange : Interchange, sim_params : SimulationParameters, **kwargs) -> Simulation:
        '''Generate an OpenMM Simulation instance using the Forces and Integrator defined for the ensemble of choice'''
        integrator = self.integrator(sim_params)
        sim = interchange.to_openmm_simulation(integrator, **kwargs)
        desc_str = f'Created {self.desc} Simulation with {integrator.__class__.__name__}'

        forces = self.forces(sim_params)
        if forces:
            for force in forces: # add forces one-by-one, per documentation in source (https://docs.openforcefield.org/projects/interchange/en/stable/_modules/openff/interchange/components/interchange.html#Interchange.to_openmm_simulation)
                sim.system.addForce(force)
            sim.context.reinitialize(preserveState=True) # reinitialize to ensure changes stick (see https://github.com/openmm/openmm/wiki/Frequently-Asked-Questions#why-does-it-ignore-changes-i-make-to-a-system-or-force) 

            force_str = ', '.join(force.__class__.__name__ for force in forces)
            desc_str = f'{desc_str} and {force_str} forces'
        LOGGER.info(desc_str)
        
        return sim

# CONCRETE IMPLEMENTATIONS
class NVESimulationFactory(EnsembleSimulationFactory):
    ensemble = 'NVE'
    ensemble_name = 'microcanonical'

    def integrator(self, sim_params: SimulationParameters) -> Integrator:
        return VerletIntegrator(sim_params.timestep)
    
    def forces(self, sim_params: SimulationParameters) -> Optional[Iterable[Force]]:
        return None
    
class NVTSimulationFactory(EnsembleSimulationFactory): # TODO : add implementation support for Andersen and Nose-Hoover thermostats (added to forces instead)
    ensemble = 'NVT'
    ensemble_name = 'canonical'

    def integrator(self, sim_params: SimulationParameters) -> Integrator:
        return LangevinMiddleIntegrator(sim_params.temperature, sim_params.friction_coeff, sim_params.timestep)
    
    def forces(self, sim_params: SimulationParameters) -> Optional[Iterable[Force]]:
        return None
    
class NPTSimulationFactory(EnsembleSimulationFactory):
    ensemble = 'NPT'
    ensemble_name = 'isothermal-isobaric'

    def integrator(self, sim_params: SimulationParameters) -> Integrator:
        return LangevinMiddleIntegrator(sim_params.temperature, sim_params.friction_coeff, sim_params.timestep)
    
    def forces(self, sim_params: SimulationParameters) -> Optional[Iterable[Force]]:
        return [MonteCarloBarostat(sim_params.pressure, sim_params.temperature, sim_params.barostat_freq)]
    
  