'''Unit and behavioral tests for thermodynamic parameter serialization'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import pytest

from typing import Optional, Type
from tempfile import NamedTemporaryFile

from openmm.unit import (
    Quantity,
    kelvin,
    atmosphere,
    picosecond,
    femtosecond,
)
from openmm import (
    Force,
    AndersenThermostat,
    MonteCarloBarostat,
    MonteCarloFlexibleBarostat,
)
from openmm import (
    Integrator,
    BrownianIntegrator,
    LangevinIntegrator,
    LangevinMiddleIntegrator,
    NoseHooverIntegrator,
    VerletIntegrator,
)

from polymerist.mdtools.openmmtools.thermo import (
    Thermostat,
    ThermostatParameters,
    Barostat,
    BarostatParameters,
    Ensemble,
    ThermoParameters,
    NPHEnsembleUnsupported,
)

# EXAMPLE PARAMETER SETS (hard-coded, since pytest STILL doesn't support parameterized fixtures >:( )
## DEVNOTE: cast as lambda to create new object on each reference, rather than passing around the same one
GET_REFERENCE_TEMPERATURE = lambda : 300*kelvin
GET_EXAMPLE_THERMOSTAT_PARAMETERS = lambda : ThermostatParameters(
    temperature=GET_REFERENCE_TEMPERATURE(),
    timescale=1/picosecond,
    thermostat=Thermostat.LANGEVIN_MIDDLE,
)
GET_EXAMPLE_BAROSTAT_PARAMETERS = lambda : BarostatParameters(
    pressure=1*atmosphere,
    temperature=GET_REFERENCE_TEMPERATURE(),
    update_frequency=25,
    barostat=Barostat.MONTE_CARLO,
)

# TESTS PROPER
@pytest.mark.parametrize(
    'thermostat,barostat,force_types_expected',
    [
        (Thermostat.ANDERSEN, Barostat.MONTE_CARLO, {AndersenThermostat, MonteCarloBarostat}),
        (Thermostat.ANDERSEN, Barostat.MONTE_CARLO_FLEXIBLE, {AndersenThermostat, MonteCarloFlexibleBarostat}),
        (Thermostat.BROWNIAN, Barostat.MONTE_CARLO, {MonteCarloBarostat}),
        (Thermostat.BROWNIAN, Barostat.MONTE_CARLO_FLEXIBLE, {MonteCarloFlexibleBarostat}),
        (Thermostat.LANGEVIN, Barostat.MONTE_CARLO, {MonteCarloBarostat}),
        (Thermostat.LANGEVIN, Barostat.MONTE_CARLO_FLEXIBLE, {MonteCarloFlexibleBarostat}),
        (Thermostat.LANGEVIN_MIDDLE, Barostat.MONTE_CARLO, {MonteCarloBarostat}),
        (Thermostat.LANGEVIN_MIDDLE, Barostat.MONTE_CARLO_FLEXIBLE, {MonteCarloFlexibleBarostat}),
        (Thermostat.NOSE_HOOVER, Barostat.MONTE_CARLO, {MonteCarloBarostat}),
        (Thermostat.NOSE_HOOVER, Barostat.MONTE_CARLO_FLEXIBLE, {MonteCarloFlexibleBarostat}),
    ]
)
def test_thermo_params_correct_forces(
        thermostat : Thermostat,
        barostat : Barostat,
        force_types_expected : set[Type[Force]],
    ) -> None:
    '''Verify that the chosen thermostat and barostat combos produce the expected Forces'''
    thermo_params = ThermoParameters(
        thermostat_params=ThermostatParameters(
            temperature=300*kelvin,
            timescale=1/picosecond,
            thermostat=thermostat,
        ),
        barostat_params=BarostatParameters(
            pressure=1*atmosphere,
            temperature=300*kelvin,
            update_frequency=25,
            barostat=barostat,
        )
    )
    force_types_actual = set(type(force) for force in thermo_params.forces())
    
    assert force_types_actual == force_types_expected, \
        f'For {thermostat} thermostat/{barostat} barostat combination' \
        f', expected forces {force_types_expected}, but got {force_types_actual} '

@pytest.mark.parametrize(
    'thermostat,barostat,integrator_type_expected',
    [
        (Thermostat.ANDERSEN, Barostat.MONTE_CARLO, VerletIntegrator),
        (Thermostat.ANDERSEN, Barostat.MONTE_CARLO_FLEXIBLE, VerletIntegrator),
        (Thermostat.BROWNIAN, Barostat.MONTE_CARLO, BrownianIntegrator),
        (Thermostat.BROWNIAN, Barostat.MONTE_CARLO_FLEXIBLE, BrownianIntegrator),
        (Thermostat.LANGEVIN, Barostat.MONTE_CARLO, LangevinIntegrator),
        (Thermostat.LANGEVIN, Barostat.MONTE_CARLO_FLEXIBLE, LangevinIntegrator),
        (Thermostat.LANGEVIN_MIDDLE, Barostat.MONTE_CARLO, LangevinMiddleIntegrator),
        (Thermostat.LANGEVIN_MIDDLE, Barostat.MONTE_CARLO_FLEXIBLE, LangevinMiddleIntegrator),
        (Thermostat.NOSE_HOOVER, Barostat.MONTE_CARLO, NoseHooverIntegrator),
        (Thermostat.NOSE_HOOVER, Barostat.MONTE_CARLO_FLEXIBLE, NoseHooverIntegrator),
    ]
)

def test_thermo_params_correct_integrator(
        thermostat : Thermostat,
        barostat : Barostat,
        integrator_type_expected : Type[Integrator],
    ) -> None:
    '''Verify that the chosen thermostat and barostat combos produce the expected Forces'''
    thermo_params = ThermoParameters(
        thermostat_params=ThermostatParameters(
            temperature=300*kelvin,
            timescale=1/picosecond,
            thermostat=thermostat,
        ),
        barostat_params=BarostatParameters(
            pressure=1*atmosphere,
            temperature=300*kelvin,
            update_frequency=25,
            barostat=barostat,
        )
    )
    integrator_type_actual = type(thermo_params.integrator(time_step=1*femtosecond))

    assert integrator_type_actual == integrator_type_expected, \
        f'For {thermostat} thermostat/{barostat} barostat combination' \
        f', expected integrator {integrator_type_expected}, but got {integrator_type_actual} '

@pytest.mark.parametrize(
    'thermostat_params,barostat_params,ensemble_expected',
    [
        (None, None, Ensemble.NVE),
        pytest.param(
            None, GET_EXAMPLE_BAROSTAT_PARAMETERS(), None,
            marks=pytest.mark.xfail(
                raises=NPHEnsembleUnsupported, 
                reason='OpenMM does not support barostat action without thermostatting',
                strict=True,
            )
        ),
        (GET_EXAMPLE_THERMOSTAT_PARAMETERS(), None, Ensemble.NVT),
        (GET_EXAMPLE_THERMOSTAT_PARAMETERS(), GET_EXAMPLE_BAROSTAT_PARAMETERS(), Ensemble.NPT),
    ]
)
def test_thermo_params_ensemble(
        thermostat_params : ThermostatParameters,
        barostat_params : BarostatParameters,
        ensemble_expected : Ensemble,
    ):
    '''Verify that the correct named thermodynamic ensemble is identified based on the provided thermostat-barostat combo'''
    thermo_params = ThermoParameters(
        thermostat_params=thermostat_params,
        barostat_params=barostat_params,
    )
    assert thermo_params.ensemble == ensemble_expected, \
        f'Expected parameter to identify {ensemble_expected}, received {thermo_params.ensemble} instead'

@pytest.mark.parametrize(
    'barostat_temperature',
    [
        None, # unset temperature
        GET_REFERENCE_TEMPERATURE() - 100*kelvin, # too low temperature
        GET_REFERENCE_TEMPERATURE() + 100*kelvin, # too high temperature
        GET_REFERENCE_TEMPERATURE(),              # just right temperature
        GET_REFERENCE_TEMPERATURE()._value, # temp w/ forgotten units (will be upconverted)
    ]
)
def test_thermo_params_temperature_coupling(barostat_temperature : Optional[Quantity]) -> None:
    '''Test that thermostat temperature is mirrored to barostat termpature, when both are provided'''
    # NOTE: termpatures must be set BEFORE initializing ThermoParamters to fatihfully test validation 
    thermostat_parameters = GET_REFERENCE_TEMPERATURE()
    thermostat_parameters.temperature = GET_REFERENCE_TEMPERATURE() # fix thermostat at constant reference
    
    barostat_parameters = GET_EXAMPLE_BAROSTAT_PARAMETERS()
    barostat_parameters.temperature = barostat_temperature # set passed temperature to test if it will be matched to the reference
    
    # reminder that by this point, temperatures should be passed as-specified, since the validation happens in __post_init__
    thermo_params = ThermoParameters(
        thermostat_params=thermostat_parameters,
        barostat_params=barostat_parameters,
    )
    
    assert thermo_params.barostat_params.temperature == thermo_params.thermostat_params.temperature, \
        f'Temperature coupling improperly enforced, ended up with thermostat temperature of' \
        f'{thermo_params.thermostat_params.temperature}, but barostat temperature of {thermo_params.barostat_params.temperature}'
    
@pytest.mark.parametrize(
    'thermostat_params,barostat_params',
    [
        (None, None),
        pytest.param(
            None, GET_EXAMPLE_BAROSTAT_PARAMETERS(),
            marks=pytest.mark.xfail(
                raises=NPHEnsembleUnsupported, 
                reason='OpenMM does not support barostat action without thermostatting',
                strict=True,
            )
        ),
        (GET_EXAMPLE_THERMOSTAT_PARAMETERS(), None),
        (GET_EXAMPLE_THERMOSTAT_PARAMETERS(), GET_EXAMPLE_BAROSTAT_PARAMETERS()),
    ]
)    

def test_thermo_params_JSON_IO(
        thermostat_params : ThermostatParameters,
        barostat_params : BarostatParameters,
    ) -> None:
    '''Test that ThermoParameters can serialize and deserialize itself and all its components to and from JSON losslessly'''
    thermo_params = ThermoParameters(
        thermostat_params=thermostat_params,
        barostat_params=barostat_params,
    )
    with NamedTemporaryFile(prefix='thermo_params', suffix='.json') as tmpfile:
        thermo_params.to_file(tmpfile.name)
        thermo_params_from_json = ThermoParameters.from_file(tmpfile.name)

    assert thermo_params == thermo_params_from_json, f'Could not perform {ThermoParameters.__class__.__name__} JSON I/O losslessly'