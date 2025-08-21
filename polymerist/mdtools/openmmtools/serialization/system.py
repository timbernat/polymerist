'''For handling serialization and description of OpenMM systems'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from pathlib import Path
from openmm import System, XmlSerializer

from ....genutils.fileutils.pathutils import allow_string_paths


@allow_string_paths
def serialize_system(sys_path : Path, system : System) -> None:
    '''For saving an existing OpenMM System to file'''
    with sys_path.open('w') as file:
        file.write(XmlSerializer.serialize(system))