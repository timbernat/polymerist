'''For converting macroscopic parameters (such as concentration, bulk density, etc) into microscopic parameters for simulations'''

from math import ceil
from rdkit.Chem import Descriptors, Mol

from openmm.unit import gram, centimeter, mole
from openmm.unit import Quantity, AVOGADRO_CONSTANT_NA

from ...unitutils.dimensions import is_volume
from ...unitutils.interop import allow_openff_units


def molecular_weight(rdmol : Mol, exact_isotopes : bool=False) -> Quantity:
    '''Compute the molecular weight of an RDKit Mol'''
    MW = Descriptors.ExactMolWt(rdmol) if exact_isotopes else Descriptors.MolWt(rdmol)
    return MW * (gram / mole) # attach units

@allow_openff_units
def number_density(density : Quantity, MW : Quantity) -> Quantity:
    '''
    Determine the number of solvent molecules per unit volume from known physical constants
    For best results, provide arguments as Quantities with associated units
    '''
    return (density / MW) * AVOGADRO_CONSTANT_NA

def num_mols_in_box(mol : Mol, box_vol : Quantity, density : Quantity, exact_isotopes : bool=False) -> int:
    '''Return the number of particles/molecules needed to fill a box of given volume to the specified density'''
    assert(is_volume(box_vol.unit))

    MW = molecular_weight(mol, exact_isotopes=exact_isotopes)
    rho_n = number_density(density=density, MW=MW)

    return ceil(box_vol * rho_n)