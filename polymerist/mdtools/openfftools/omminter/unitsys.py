'''For handling interconversion between the OpenMM and OpenFF (Pint) unit engines'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Callable, TypeVar

ReturnType = TypeVar('ReturnType') 

from openmm.unit import Quantity as OpenMMQuantity
from pint import Quantity as PintQuantity # this is also the base class for all OpenFF-style units

from openff.units.openmm import (
    from_openmm as openmm_to_openff,
    to_openmm as openff_to_openmm,
)


def allow_openmm_units(funct : Callable[..., ReturnType]) -> Callable[..., ReturnType]:
    '''Allow a Callable which expects any of its args to be OpenFF Quantities to also accept equivalent OpenMM Quantites'''
    def wrapper(*args, **kwargs) -> ReturnType:
        new_args = [
            openmm_to_openff(arg) if isinstance(arg, OpenMMQuantity) else arg
                for arg in args
        ]

        new_kwargs = {
            key : openmm_to_openff(kwarg) if isinstance(kwarg, OpenMMQuantity) else kwarg
                for key, kwarg in kwargs.items()
        }

        return funct(*new_args, **new_kwargs)
    return wrapper

def allow_openff_units(funct : Callable[..., ReturnType]) -> Callable[..., ReturnType]:
    '''Allow a Callable which expects any of its args to be OpenMM Quantities to also accept equivalent OpenFF Quantites'''
    def wrapper(*args, **kwargs) -> ReturnType:
        new_args = [
            openff_to_openmm(arg) if isinstance(arg, PintQuantity) else arg
                for arg in args
        ]

        new_kwargs = {
            key : openff_to_openmm(kwarg) if isinstance(kwarg, PintQuantity) else kwarg
                for key, kwarg in kwargs.items()
        }

        return funct(*new_args, **new_kwargs)
    return wrapper

# TODO : add decorators which allows OUTPUTS to have flexible units