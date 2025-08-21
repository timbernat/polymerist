'''Utilities for extending the Python-JSON interface'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from .jsonify import JSONifiable, make_jsonifiable
from .serialize import (
    JSONSerializable,
    TypeSerializer,
    MultiTypeSerializer,
    PathSerializer,
    QuantitySerializer
)
from .update import append_to_json, JSONDict
