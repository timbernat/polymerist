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
)
ureg = UnitRegistry()

from polymerist.unitutils.dimensions import (
    hasunits,
    strip_units,
    is_volume,
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