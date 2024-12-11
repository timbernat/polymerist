'''For converting macroscopic parameters (such as concentration, bulk density, etc) into microscopic parameters for simulations'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Union

from math import ceil
from rdkit.Chem import Descriptors, Mol

from openmm.unit import gram, centimeter, mole
from openmm.unit import Quantity, AVOGADRO_CONSTANT_NA

from openff.toolkit import Molecule, Topology
from openff.units import Quantity as OFFQuantity

from ....unitutils.dimensions import is_volume
from ..unitsys import allow_openff_units, openff_to_openmm


# MASS
def molecular_weight(rdmol : Mol, exact_isotopes : bool=False) -> Quantity:
    '''Compute the molecular weight of an RDKit Mol'''
    MW = Descriptors.ExactMolWt(rdmol) if exact_isotopes else Descriptors.MolWt(rdmol)
    return MW * (gram / mole) # attach units

def offmol_mass(offmol : Molecule, as_openmm : bool=False) -> Union[Quantity, OFFQuantity]:
    '''Compute the molecular mass of an OpenFF Molecule (in Daltons)
    Returns mass with OpenFF-style units, or with OpenMM-style units if as_openmm=True'''
    mass = sum(atom.mass for atom in offmol.atoms)
    if as_openmm:
        return openff_to_openmm(mass)
    return mass

def offtop_mass(offtop : Topology, as_openmm : bool=False) -> Union[Quantity, OFFQuantity]:
    '''Compute the molecular mass of an OpenFF Topology (in Daltons)
    Returns mass with OpenFF-style units, or with OpenMM-style units if as_openmm=True'''
    mass = sum(offmol_mass(offmol, as_openmm=False) for offmol in offtop.molecules)
    if as_openmm:
        return openff_to_openmm(mass) # perform conversion at end, to avoid multiple intermediate conversions and avoid initializing sum() with Quantity
    return mass

# DENSITY
@allow_openff_units
def number_density(density : Quantity, MW : Quantity) -> Quantity:
    '''
    Determine the number of solvent molecules per unit volume from known physical constants
    For best results, provide arguments as Quantities with associated units
    '''
    return (density / MW) * AVOGADRO_CONSTANT_NA

# NUMBER
@allow_openff_units
def num_mols_in_box(mol : Union[Mol, Molecule, Topology], box_vol : Quantity, density : Quantity) -> int:
    '''Return the number of particles/molecules needed to fill a box of given volume to the specified density'''
    assert(is_volume(box_vol.unit))

    if isinstance(mol, Mol):
        MW = molecular_weight(mol, exact_isotopes=False)
    elif isinstance(mol, Molecule):
        MW = offmol_mass(mol, as_openmm=True)
    elif isinstance(mol, Topology):
        MW = offtop_mass(mol, as_openmm=True)
    else:
        raise TypeError(f'No molecular weight calculations implemented for object of type "{type(mol).__name__}"')

    return ceil(box_vol * number_density(density=density, MW=MW))