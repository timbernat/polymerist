'''Tools for editing and manipulating images, and image colors, sizes, and representations'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Union

import numpy as np
from io import BytesIO

import PIL
from PIL.Image import Image


# CONVERSION FUNCTIONS
def img_from_bytes(img_bytes : bytearray) -> Image:
    '''Load an image from a bytestream'''
    return PIL.Image.open(BytesIO(img_bytes))

def img_to_array(image : Image, encoding : str='RGB') -> np.ndarray:
    '''Convert an image to a numpy array'''
    return np.asarray(image.convert(encoding))

# BOUNDING BOXES AND BACKGROUND REMOVAL
def get_axis_bounds(image : Image, bg_color : Union[int, tuple[int]]) -> tuple[tuple[int, int], tuple[int, int]]:
    '''Takes an image and returns a pair of tuples containing the min and max non-background
    coordinates along the x- and y-axis (respectively) for some chosen background color'''
    y, x, ch = np.where((np.asarray(image) != bg_color)) # TODO : generalize this to work for greyscale (too many dimensions)
    x_bounds = (x.min(), x.max())
    y_bounds = (y.min(), y.max())

    return x_bounds, y_bounds

def get_tight_bbox(image : Image, bg_color : Union[int, tuple[int]]) -> tuple[int, ...]:
    '''Takes an image and returns a 4-tuple of the coordinates of the tight bounding box based on a choice of background color'''
    (x_min, x_max), (y_min, y_max) = get_axis_bounds(image, bg_color=bg_color)
    tight_bbox = (x_min, y_min, x_max, y_max)

    return tight_bbox

def crop_borders(image : Image, bg_color : Union[int, tuple[int]]) -> Image:
    '''Takes an image and returns an image with all extraneous borders of a particular background color cropped off'''
    return image.crop(get_tight_bbox(image, bg_color=bg_color)) 