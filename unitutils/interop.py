'''Decorators for handling interconversion between the OpenMM and OpenFF (Pint) unit engines'''

from typing import Callable, TypeVar
R = TypeVar('R')   # for representing generic return values
Q = TypeVar('Q')   # for representing generic Quantity-like objects

from pint import Quantity as PintQuantity # this is also the base class for all OpenFF-style units
from openmm.unit import Quantity
from openff.units.openmm import (
    from_openmm as openmm_to_openff,
    to_openmm as openff_to_openmm,
)


def allow_openmm_units(funct : Callable[[Q], R]) -> Callable[[Q], R]:
    '''Allow a Callable which expects ALL of its args to be OpenFF Quanitities to also accept equivalent OpenMM Quanitites'''
    def wrapper(*args, **kwargs) -> R:
        new_args = [
            openmm_to_openff(arg) if isinstance(arg, Quantity) else arg
                for arg in args
        ]

        new_kwargs = {
            key : openmm_to_openff(kwarg) if isinstance(kwarg, Quantity) else kwarg
                for key, kwarg in kwargs.items()
        }

        return funct(*new_args, **new_kwargs)
    return wrapper

def allow_openff_units(funct : Callable[[Q], R]) -> Callable[[Q], R]:
    '''Allow a Callable which expects ALL of its args to be OpenMM Quanitities to also accept equivalent OpenFF Quanitites'''
    def wrapper(*args, **kwargs) -> R:
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

