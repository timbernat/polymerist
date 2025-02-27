'''Representations and conversion methods for colors'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import TypeAlias, Union

RGB : TypeAlias = tuple[int, int, int]
RGBA : TypeAlias = tuple[int, int, int, int]
Color = Union[hex, RGB, RGBA] # TODO : add typecheck method

NAMED_COLORS : dict[str, RGB] = {
    # greyscale
    'BLACK' : (0, 0, 0),
    'WHITE' : (255, 255, 255),
    # primary
    'RED'   : (255, 0, 0),
    'GREEN' : (0, 255, 0),
    'BLUE'  : (0, 0, 255),
}
# dynamically register named colors at the module level
for color_name, color_code in NAMED_COLORS.items():
    globals()[color_name] = color_code