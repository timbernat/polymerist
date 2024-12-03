'''For dynamically determining and cataloging which SMIRNOFF-copatible force fields are installed (and accompanying functionality) are available'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Optional
from pathlib import Path

from ...genutils.importutils.dependencies import modules_installed


# Force field and ToolkitWrapper reference
FFDIR : Optional[Path] = None
if modules_installed('openff.toolkit'):
    from openforcefields import get_forcefield_dirs_paths
    
    FFDIR = Path(get_forcefield_dirs_paths()[0]) # Locate path where OpenFF forcefields are installed

FF_DIR_REGISTRY  : dict[Path, Path] = {}
FF_PATH_REGISTRY : dict[Path, Path] = {}
if modules_installed('openforcefields'):
    from openff.toolkit.typing.engines.smirnoff.forcefield import _get_installed_offxml_dir_paths
    
    for ffdir_str in _get_installed_offxml_dir_paths():
        ffdir = Path(ffdir_str)
        ffdir_name = ffdir.parent.stem

        FF_DIR_REGISTRY[ ffdir_name]  = ffdir
        FF_PATH_REGISTRY[ffdir_name] = [path for path in ffdir.glob('*.offxml')]