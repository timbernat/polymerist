'''API for selecting thermostat and barostat actions which realize particular thermodynamic ensembles'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import logging
LOGGER = logging.getLogger(__name__)

from typing import Iterable, Optional, Union
from dataclasses import dataclass, field
from enum import Enum, StrEnum

from openmm.unit import Quantity, picosecond, kelvin, atmosphere
from openmm import Force, Integrator
from openmm import (
    LangevinMiddleIntegrator,
    LangevinIntegrator,
    BrownianIntegrator,
    NoseHooverIntegrator,
    AndersenThermostat,
    VerletIntegrator,
)
from openmm import (
    MonteCarloBarostat,
    MonteCarloFlexibleBarostat,
)

from polymerist.genutils.fileutils.jsonio.jsonify import make_jsonifiable
from polymerist.genutils.fileutils.jsonio.serialize import (
    MultiTypeSerializer,
    QuantitySerializer,
    enum_serializer_factory,
)


# THERMOSTATS
class Thermostat(Enum):
    '''Common thermostats which OpenMM implements and which adhere to the interface defined here'''
    ANDERSEN = AndersenThermostat
    BROWNIAN = BrownianIntegrator
    LANGEVIN = LangevinIntegrator
    LANGEVIN_MIDDLE = LangevinMiddleIntegrator
    NOSE_HOOVER = NoseHooverIntegrator
    
    # aliases for convenience
    LANGEVINMIDDLE = LangevinMiddleIntegrator
    NOSEHOOVER = NoseHooverIntegrator
ThermostatSerializer = enum_serializer_factory(Thermostat)

@make_jsonifiable(type_serializer=MultiTypeSerializer(QuantitySerializer, ThermostatSerializer))
@dataclass
class ThermostatParameters:
    '''Interface for initializing a constant-temperature simulation'''
    temperature : Quantity
    timescale : Quantity = field(default_factory=lambda : 1 * picosecond**-1)
    thermostat : Union[str, Thermostat] = Thermostat.LANGEVIN_MIDDLE

    def __post_init__(self):
        if isinstance(self.thermostat, str): # DEVNOTE: allowing strings is forgiving to users oblivious to enums here
            self.thermostat = Thermostat[self.thermostat.upper()]

    def forces(self) -> Iterable[Force]:
        '''The forces required to realized the desired thermostat'''
        if self.thermostat == Thermostat.ANDERSEN: # DEVNOTE: to my unending frustration, Andersen is the ONLY thermostat implemented as a Force (not an integrator)!
            return (Thermostat.ANDERSEN.value(self.temperature, self.timescale),)
        return tuple()

    def integrator(self, time_step : Quantity) -> Integrator:
        '''The integrator required to realized the desired thermostat'''
        if self.thermostat == Thermostat.ANDERSEN: # DEVNOTE: to my unending frustration, Andersen is the ONLY thermostat implemented as a Force (not an integrator)!
            return VerletIntegrator(time_step)
        return self.thermostat.value(self.temperature, self.timescale, time_step)

# BAROSTATS
class Barostat(Enum):
    '''Common barostats which OpenMM implements and which adhere to the interface defined here'''
    # DEVNOTE: don't support Anisotropic or Membrane Barostats, since their call signatures don't fit the same interface as the 
    # "vanilla" MonteCarloBarostat (vectored pressure and scale toggles for the former, non-optional surface tension for the latter)
    MONTE_CARLO = MonteCarloBarostat
    MONTE_CARLO_FLEXIBLE = MonteCarloFlexibleBarostat

    # aliases for convenience
    MC = MonteCarloBarostat
    MONTECARLO = MonteCarloBarostat
    FLEXIBLE = MonteCarloFlexibleBarostat
BarostatSerializer = enum_serializer_factory(Barostat)

@make_jsonifiable(type_serializer=MultiTypeSerializer(QuantitySerializer, BarostatSerializer))
@dataclass
class BarostatParameters:
    '''Interface for initializing a constant-pressure simulation'''
    pressure : Quantity
    temperature : Optional[Quantity] = None # DEVNOTE: made NoneType to allow "lazy" passing when coupled with thermostat 
    update_frequency : int = 25
    barostat : Union[str, Barostat] = Barostat.MONTE_CARLO

    def __post_init__(self):
        if isinstance(self.barostat, str): # DEVNOTE: allowing strings is forgiving to users oblivious to enums here
            self.barostat = Barostat[self.barostat.upper()]

    def forces(self) -> Iterable[Force]:
        '''The forces required to realized the desired barostat'''
        if self.temperature is None:
            raise ValueError('Barostat coupling temperature unset')
        return (self.barostat.value(self.pressure, self.temperature, self.update_frequency),)

    def integrator(self, time_step : Quantity) -> Integrator:
        '''The integrator required to realized the desired barostat'''
        # DEVNOTE: just here to keep interface constant; since the NPH ensemble is not supported,
        # there is never a need for the barostat to provide an integrator
        return tuple()
    
# ENSEMBLES
class NPHEnsembleUnsupported(ValueError):
    '''
    Raised when a user attempts to initialize ThermoParameters with a barostat but not thermostat
    Would cause OpenMM to produce incorrect results (https://docs.openmm.org/latest/userguide/application/02_running_sims.html#pressure-coupling)
    '''
    def __init__(self, msg : str=f'NPH ensemble not supported; either add a thermostat or remove a barostat from thermodynamic parameters', *args, **kwargs):
        super().__init__(msg, *args, **kwargs)

class Ensemble(StrEnum):
    '''Common thermodynamic ensembles which are realizable with the interfaces provided here'''
    NVE = 'microcanonical'
    NVT = 'canonical'
    NPT = 'isothermal-isobaric'
    
@make_jsonifiable
@dataclass
class ThermoParameters:
    '''Encapsulation for initializing the OpenMM forces and integrator which realize a particular thermodynamic ensemble'''
    thermostat_params : Optional[ThermostatParameters] = None
    barostat_params   : Optional[BarostatParameters  ] = None
    
    def __post_init__(self) -> None:
        if self.barostat_params is not None:
            if self.thermostat_params is None:
                raise NPHEnsembleUnsupported
            
            if self.barostat_params.temperature != self.thermostat_params.temperature:
                LOGGER.warning(f'Adjusting Barostat temperature from "{self.barostat_params.temperature}" to {self.thermostat_params.temperature} to maintain temperature coupling w/ thermostat')
                self.barostat_params.temperature = self.thermostat_params.temperature

    @property
    def ensemble(self) -> Ensemble:
        '''The standard name of the thermodynamic ensemble being implemented here'''
        if self.thermostat_params is None:
            if self.barostat_params is not None:
                raise NPHEnsembleUnsupported
            return Ensemble.NVE
        else:
            if self.barostat_params is not None:
                return Ensemble.NPT
            return Ensemble.NVT
        
    def describe_ensemble(self) -> str:
        '''Verbal description of ensemble'''
        return f'{self.ensemble.name} ({self.ensemble.value.capitalize()}) ensemble'
        
    def integrator(self, time_step : Quantity) -> Integrator:
        '''Specify how to integrate forces in each timestep'''
        if self.thermostat_params:
            integrator = self.thermostat_params.integrator(time_step)
        else:
            integrator = VerletIntegrator(time_step)

        LOGGER.info(f'Created {integrator.__class__.__name__} for {self.describe_ensemble()}')
        
        return integrator

    def forces(self) -> Optional[Iterable[Force]]:
        '''Specify any additional force contributions to position/velocity updates'''
        forces : list[Force] = []
        if self.thermostat_params:
            forces.extend(self.thermostat_params.forces())
        if self.barostat_params:
            forces.extend(self.barostat_params.forces())

        for force in forces:
            LOGGER.info(f'Created {force.__class__.__name__} Force for {self.describe_ensemble()}')

        return forces