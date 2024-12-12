'''For curating pre-defined solvent molecules'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from pathlib import Path # TODO : reimplement "properly" "using importlib_resources (more complicated than it's worth for now)
_MODULE_PATH = Path(__path__[0])

from openff.toolkit import Molecule
from openff.units import unit as offunit

from ... import topology


def generate_water_TIP3P() -> Molecule:
    '''Helper method for creating a new TIP3P water representation from scratch'''
    TIP3P_ATOM_CHARGES = { # NOTE : units deliberately omitted here (become applied to entire charge array)
        'H' :  0.417,
        'O' : -0.843
    }

    water = Molecule.from_smiles('O')
    water.name = 'water_TIP3P'
    water.partial_charges = [TIP3P_ATOM_CHARGES[atom.symbol] for atom in water.atoms]*offunit.elementary_charge
    for atom in water.atoms:
        atom.metadata['residue_name'] = 'HOH'

    return water


# predefine water file, if not already present
_water_path = _MODULE_PATH / 'water_TIP3P.sdf'
if not _water_path.exists():
    water = generate_water_TIP3P()
    topology.package_atom_metadata(water, in_place=True)
    topology.save_molecule(_water_path, water)

# register Molcules for all registered solvents
for path in _MODULE_PATH.iterdir():
    if path.suffix == '.sdf':
        offmol = Molecule.from_file(path) # load molecules from file
        globals()[path.stem] = topology.unpackage_atom_metadata(offmol, in_place=False) # recover metadata and register solvent Molecules locally
