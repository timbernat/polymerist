'''Testing that dimensionality checking behaves as expected for both OpenMm and Pint-style unit systems'''

from typing import Any, Union
from dataclasses import dataclass

import pytest

from openmm.unit import (
    Unit as OpenMMUnit,
    Quantity as OpenMMQuantity,
    centimeter,
    second,
)
from pint import (
    Unit as PintUnit,
    Quantity as PintQuantity, # this is also the base class for all OpenFF-style units
    UnitRegistry,
    DimensionalityError,
)
ureg = UnitRegistry()
Quantity = Union[PintQuantity, OpenMMQuantity]

from polymerist.unitutils.dimensions import (
    hasunits,
    strip_units,
    is_volume,
    quantities_approx_equal
)


# Defining test cases and expected outputs
@dataclass
class UnitExample:
    '''Internal encapsulation class for indicating expected
    properties of unit-like objects in unit tests (no pun intended)'''
    value : Any
    has_units : bool
    is_a_volume : bool # changed name slightly to obviate clash with is_volume() in namespace
    
test_cases : list[UnitExample] = [
    # non-units
    UnitExample(
        value=42,
        has_units=False,
        is_a_volume=False,
    ),
    UnitExample(
        value=3.1415,
        has_units=False,
        is_a_volume=False,
    ),
    UnitExample(
        value={1,2,3},
        has_units=False,
        is_a_volume=False,
    ),
    # pure units
    UnitExample(
        value=second,
        has_units=False,
        is_a_volume=False,
    ),
    UnitExample(
        value=centimeter**3,
        has_units=False,
        is_a_volume=True, # despite being a pure unit, this should still count as a volume
    ),
    UnitExample(
        value=ureg.second,
        has_units=False,
        is_a_volume=False,
    ),
    UnitExample(
        value=ureg.foot**3,
        has_units=False,
        is_a_volume=True,
    ),
    # simple quantities
    UnitExample(
        value=1.1*second,
        has_units=True,
        is_a_volume=False,
    ),
    UnitExample(
        value=1.2*centimeter,
        has_units=True,
        is_a_volume=False,
    ),
    UnitExample(
        value=1.3*centimeter**3,
        has_units=True,
        is_a_volume=True,
    ),
    UnitExample(
        value=2.1*ureg.second,
        has_units=True,
        is_a_volume=False,
    ),
    UnitExample(
        value=2.2*ureg.centimeter,
        has_units=True,
        is_a_volume=False,
    ),
    UnitExample(
        value=2.3*ureg.centimeter**3,
        has_units=True,
        is_a_volume=True,
    ),
    # mixed quantities
    UnitExample(
        value=9.8*centimeter*second**-1,
        has_units=True,
        is_a_volume=False,
    ),
    UnitExample(
        value=9.81*ureg.centimeter*ureg.second**-1,
        has_units=True,
        is_a_volume=False,
    )
]

# Unit tests 
@pytest.mark.parametrize('unitlike, expected_output', [
    (unit_example.value, unit_example.has_units)
        for unit_example in test_cases
    ]
)
def test_hasunits(unitlike : Any, expected_output : bool) -> None:
    '''Test that objects with (and without) units are correctly identified'''
    assert hasunits(unitlike) == expected_output

@pytest.mark.parametrize('unitlike, expected_output', [
    (unit_example.value, unit_example.is_a_volume)
        for unit_example in test_cases
    ]
)
def test_hasunits(unitlike : Any, expected_output : bool) -> None:
    '''Test that objects which can (and can't) be interpreted as volumes are correctly identified as such'''
    assert is_volume(unitlike) == expected_output
    
SAMPLE_COORDS : tuple[float] = (1.23, 4.56, 7.89) # these numbers are arbitrary, but need to be consistent across tests
@pytest.mark.parametrize('coordlike', [SAMPLE_COORDS, SAMPLE_COORDS*centimeter, SAMPLE_COORDS*ureg.centimeter])
def test_strip_units(coordlike : Union[tuple, PintQuantity, OpenMMQuantity]) -> None:
    '''Test that removing units works for Pint, OpenMM, and unit-free objects'''
    assert tuple(strip_units(coordlike)) == SAMPLE_COORDS # need to re-tuplify to counteract numpy auto-conversion by pint

@pytest.mark.parametrize(
    'q1, q2, rel_tol, expected_are_equal', [
        # test nominal OpenMM inputs
        (1.0*centimeter, (1.0 + 1E-6)*centimeter, 1E-4, True),
        (1.0*centimeter, (1.0 + 1E-6)*centimeter, 1E-8, False),
        ## check within float precision
        (0.1*(centimeter/second) + 0.2*(centimeter/second), 0.3*(centimeter/second), 1E-8, True), 
        (0.1*(centimeter/second) + 0.2*(centimeter/second), 0.3*(centimeter/second), 1E-18, False),
        # test nominal Pint inputs
        (1.0*ureg.centimeter, (1.0 + 1E-6)*ureg.centimeter, 1E-4, True),
        (1.0*ureg.centimeter, (1.0 + 1E-6)*ureg.centimeter, 1E-8, False),
        ## check within float precision
        (0.1*ureg.centimeter_per_second + 0.2*ureg.centimeter_per_second, 0.3*ureg.centimeter_per_second, 1E-8, True), 
        (0.1*ureg.centimeter_per_second + 0.2*ureg.centimeter_per_second, 0.3*ureg.centimeter_per_second, 1E-18, False),
        # Pathological examples
        pytest.param(
            1.0*centimeter, (1.0 + 1E-6)*ureg.centimeter, 1E-4, True, 
            marks=pytest.mark.xfail(
                raises=TypeError, 
                reason='OpenMM and Pint units should not be directly comparable',
                strict=True,
            )
        ),
        pytest.param( # result is not symmetric in order of args, so check swapped args for good measure
            (1.0 + 1E-6)*ureg.centimeter, 1.0*centimeter, 1E-4, True, 
            marks=pytest.mark.xfail(
                raises=TypeError, 
                reason='OpenMM and Pint units should not be directly comparable',
                strict=True,
            )
        ),
        pytest.param(
            3.14, 3.1415, 1E-4, True,
            marks=pytest.mark.xfail(
                raises=TypeError,
                reason='Non-quantity objects should not be comparable',
                strict=True,
            )
        ),
        pytest.param(
            1.0*centimeter, 1.0*second, 1E-4, True,
            marks=pytest.mark.xfail(
                raises=TypeError, # Exception raised will be specific to OpenMM...
                reason='Incompatible units should not be comparable',
                strict=True,
            )
        ),
        pytest.param(
            1.0*ureg.centimeter, 1.0*ureg.second, 1E-4, True,
            marks=pytest.mark.xfail(
                raises=DimensionalityError, # ...vs Pint
                reason='Incompatible units should not be comparable',
                strict=True,
            )
        ),
    ]
)
def test_quantities_approx_equal(q1 : Quantity, q2 : Quantity, rel_tol : float, expected_are_equal : bool) -> None:
    '''Test that two quantities are approximately equal within a given tolerance'''
    assert quantities_approx_equal(q1, q2, rel_tol) == expected_are_equal