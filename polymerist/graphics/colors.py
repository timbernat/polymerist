'''Representations and conersion methods for colors'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import TypeAlias, Union

RGB : TypeAlias = tuple[int, int, int]
RGBA : TypeAlias = tuple[int, int, int, int]
Color = Union[hex, RGB, RGBA] # TODO : add typecheck method