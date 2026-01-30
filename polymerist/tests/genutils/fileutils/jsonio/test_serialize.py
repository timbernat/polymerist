'''Unit tests for auto-JSONification of dataclasses, as enabled by TypeSerializers'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import pytest

from typing import Any, Type, Union
from dataclasses import dataclass, asdict

import json
from pathlib import Path
from tempfile import NamedTemporaryFile

import numpy as np
from numpy.testing import assert_equal
from openmm.unit import Quantity, nanometer, picosecond

from polymerist.genutils.fileutils.jsonio.jsonify import make_jsonifiable
from polymerist.genutils.fileutils.jsonio.serialize import (
    TypeSerializer,
    MultiTypeSerializer,
    SetSerializer,
    PathSerializer,
    QuantitySerializer,
    NDArraySerializer,
) 
PanTypeSerializer = MultiTypeSerializer(
    SetSerializer,
    PathSerializer,
    QuantitySerializer,
    NDArraySerializer,
)

# Testing JSON encode/decode for particular types    
@pytest.mark.parametrize(
    'py_obj, ts_type',
    [
        ({'Russell', 'Whitehead', 'anything but this set'}, SetSerializer),
        (Path('temp/documents/info.txt'), PathSerializer),
        (np.arange(5), NDArraySerializer),
        (Quantity(1, nanometer), QuantitySerializer),
    ]
)
def test_type_encode(
    py_obj : Any,
    ts_type : Type[TypeSerializer],
) -> None:
    '''Test that a given type serializer can JSON-encode an input of a given type'''
    try:
        json_str : str = json.dumps(ts_type.encode(py_obj))
    except TypeError:
        pytest.fail(f'Could not encode term {py_obj!r} of type {type(py_obj)} to JSON')

@pytest.mark.parametrize(
    'json_obj, ts_type',
    [
        ("{'Russell', 'Whitehead', 'anything but this set'}", SetSerializer),
        ('temp/documents/info.txt', PathSerializer),
        ({'array': [0, 1, 2, 3, 4], 'dtype': 'int64'}, NDArraySerializer),
        pytest.param(
            [0, 1, 2, 3, 4], NDArraySerializer,
            marks=pytest.mark.xfail(
                raises=TypeError,
                reason='Array decode requires dict encoding array and dtype',
                strict=True,
            )
        ),
        ({'value': 1, 'unit': 'nanometer'}, QuantitySerializer),
        pytest.param(
            '1 nanometer', QuantitySerializer,
            marks=pytest.mark.xfail(
                raises=TypeError,
                reason='Quantity decode requires dict encoding magnitude and unit',
                strict=True,
            )
        ),
    ]
)
def test_type_decode(
    json_obj : Any,
    ts_type : Type[TypeSerializer],
) -> None:
    '''Test that a given type serializer can JSON-decode an input of a given type'''
    decoded_obj = ts_type.decode(json_obj)
    assert issubclass(type(decoded_obj), ts_type.python_type)


# Testing fidelity of dataclass JOSNification
## Set
@make_jsonifiable(type_serializer=SetSerializer)
@dataclass
class DummySetContainer:
    single_field : int
    list_field : list[int]
    dict_field : dict[str, Union[str, set]]

@make_jsonifiable
@dataclass
class DummySetContainerIncomplete:
    single_field : int
    list_field : list[int]
    dict_field : dict[str, Union[str, set]]

## Path
@make_jsonifiable(type_serializer=PathSerializer)
@dataclass
class DummyPathContainer:
    single_field : int
    list_field : list[int]
    dict_field : dict[str, Union[str, Path]]

@make_jsonifiable
@dataclass
class DummyPathContainerIncomplete:
    single_field : int
    list_field : list[int]
    dict_field : dict[str, Union[str, Path]]

## Numpy array
@make_jsonifiable(type_serializer=NDArraySerializer)
@dataclass
class DummyNDArrayContainer:
    single_field : np.ndarray
    list_field : list[np.ndarray]
    dict_field : dict[str, Union[str, np.ndarray]]

@make_jsonifiable
@dataclass
class DummyNDArrayContainerIncomplete:
    single_field : np.ndarray
    list_field : list[np.ndarray]
    dict_field : dict[str, Union[str, np.ndarray]]

## Quantity
@make_jsonifiable(type_serializer=QuantitySerializer)
@dataclass
class DummyQuantityContainer:
    single_field : Quantity
    list_field : list[Quantity]
    dict_field : dict[str, Union[str, Quantity]]

@make_jsonifiable
@dataclass
class DummyQuantityContainerIncomplete:
    single_field : Quantity
    list_field : list[Quantity]
    dict_field : dict[str, Union[str, Quantity]]
    
## Multiple types at once
@make_jsonifiable(type_serializer=PanTypeSerializer)
@dataclass
class DummyMultiTypeContainer:
    name : str
    aliases : set[str]
    idx : int
    file : Path
    times : np.ndarray
    length : Quantity

@make_jsonifiable
@dataclass
class DummyMultiTypeContainerIncomplete:
    name : str
    aliases : set[str]
    idx : int
    file : Path
    times : np.ndarray
    length : Quantity
    
## Tests proper
def jsonifiable_dataclass_test_inputs() -> tuple[tuple[Type, dict[str, Any]]]:
    '''Compactly produce example test inputs for JSONifiable dataclasses'''
    DUMMY_ARGS : dict[type, dict[str, Any]] = {
        set : {
            'single_field' : {1,2,3},
            'list_field' : [{1,2,3}, set('abc')],
            'dict_field' : {
                'numbs' : {4,5,6},
                'letts' : set('def'),
            },
        },
        Path : {
            'single_field' : Path('docs/foo.txt'),
            'list_field' : [Path('docs/bar.txt'), Path('docs/bar.txt')],
            'dict_field' : {
                'foo' : Path('docs/foo.txt'),
                'bar' : Path('docs/bar.txt'),
            },
        },
        np.ndarray : {
            'single_field' : np.linspace(0, 1, 4),
            'list_field' : [np.array([0,1,2]), np.arange(5)],
            'dict_field' : {
                'arange' : np.arange(5),
                'linspace' : np.linspace(0, 1, 4),
            },
        },
        Quantity : {
            'single_field' : 420.69*nanometer,
            'list_field' : [3.14*picosecond, 0.5772*nanometer],
            'dict_field' : {
                'length' : 1*nanometer,
                'duration' : 5*picosecond,
            },
        },
        Any : {
            'name': 'Dracula',
            'aliases': {'D', 'Drac', 'Dracul'},
            'idx': 42,
            'file': Path('temp/documents/info.txt'),
            'times': np.array([0, 1, 2, 3, 4]),
            'length': Quantity(1*picosecond),
        }
    }
    
    return (
        (DummySetContainer, DUMMY_ARGS[set]),
        pytest.param(
            DummySetContainerIncomplete, DUMMY_ARGS[set],
            marks=pytest.mark.xfail(
                raises=TypeError,
                reason='Forgot to include TypeSerializer in dataclass definition',
                strict=True,
            )
        ),
        (DummyPathContainer, DUMMY_ARGS[Path]),
        pytest.param(
            DummyPathContainerIncomplete, DUMMY_ARGS[Path],
            marks=pytest.mark.xfail(
                raises=TypeError,
                reason='Forgot to include TypeSerializer in dataclass definition',
                strict=True,
            )
        ),
        (DummyNDArrayContainer, DUMMY_ARGS[np.ndarray]),
        pytest.param(
            DummyNDArrayContainerIncomplete, DUMMY_ARGS[np.ndarray],
            marks=pytest.mark.xfail(
                raises=TypeError,
                reason='Forgot to include TypeSerializer in dataclass definition',
                strict=True,
            )
        ),
        (DummyQuantityContainer, DUMMY_ARGS[Quantity]),
        pytest.param(
            DummyQuantityContainerIncomplete, DUMMY_ARGS[Quantity],
            marks=pytest.mark.xfail(
                raises=TypeError,
                reason='Forgot to include TypeSerializer in dataclass definition',
                strict=True,
            )
        ),
        (DummyMultiTypeContainer, DUMMY_ARGS[Any]),
        pytest.param(
            DummyMultiTypeContainerIncomplete, DUMMY_ARGS[Any],
            marks=pytest.mark.xfail(
                raises=TypeError,
                reason='Forgot to include TypeSerializer in dataclass definition',
                strict=True,
            )
        ),
    )
    
@pytest.mark.parametrize(
    'container_type, args',
    jsonifiable_dataclass_test_inputs()
)
def test_jsonify_dataclass_serialize(
    container_type : Type,
    args : dict[str, Any],
) -> None:
    '''
    Test that custom TypeSerializers passed to dataclasses
    wrapped with `make_jsonifiable` are able to actually 
    serialize the target types to JSON
    '''
    test_obj = container_type(**args)
    with NamedTemporaryFile('r', suffix='.json') as file:
        test_obj.to_file(file.name) # DEV: no assert - just checking that a TypeError is not raised here

@pytest.mark.parametrize(
    'container_type, args',
    jsonifiable_dataclass_test_inputs()
)
def test_jsonify_dataclass_deserialize(
    container_type : Type,
    args : dict[str, Any],
) -> None:
    test_obj = container_type(**args)
    with NamedTemporaryFile('r', suffix='.json') as file:
        test_obj.to_file(file.name)
        test_obj_loaded = container_type.from_file(file.name)
        
        assert(asdict(test_obj), asdict(test_obj_loaded))