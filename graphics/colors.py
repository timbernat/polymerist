'''Representations and conersion methods for colors'''

from typing import TypeAlias, Union


RGB : TypeAlias = tuple[int, int, int]
RGBA : TypeAlias = tuple[int, int, int, int]
Color = Union(hex, RGB, RGBA)