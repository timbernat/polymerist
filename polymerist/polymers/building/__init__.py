'''
Tools for building polymer conformers out of monomer SMARTS fragments
Currently restricted to building linear homopolymers and periodic block copolymers
'''

from ...genutils.importutils.dependencies import modules_installed, MissingPrerequisitePackage

if not modules_installed('mbuild'):
    raise MissingPrerequisitePackage(
        importing_package_name=__spec__.name,
        use_case='Polymer building',
        install_link='https://mbuild.mosdef.org/en/stable/getting_started/installation/installation.html',
        dependency_name='mbuild',
        dependency_name_formal='mBuild',
    )
    
from .linear import build_linear_polymer
from .mbconvert import (
    mbmol_from_mono_rdmol, mbmol_to_rdmol,
    mbmol_to_openmm_pdb, mbmol_to_rdkit_pdb,
)