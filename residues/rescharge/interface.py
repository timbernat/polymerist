'''Interfaces between residue-charge calculation methods and OpenFF (or other external) tools'''

from typing import ClassVar
from dataclasses import dataclass

from openff.toolkit import Molecule

from .rctypes import ChargesByResidue
from .calculation import apply_residue_charges
from ...openfftools.pcharge import MolCharger
from ...genutils.decorators.functional import optional_in_place


@dataclass
class LibraryCharger(MolCharger):
    '''Charger class for applying library charges onto residue-mapped Molecules'''
    charges_by_residue : ChargesByResidue
    CHARGING_METHOD : ClassVar[str] = 'RCT_residue_charges'

    @optional_in_place
    def _charge_molecule(self, uncharged_mol : Molecule) -> None:
        apply_residue_charges(uncharged_mol, self.charges_by_residue, in_place=True) # TOSELF: class is defined at runtime, so the fact that this functions is defined later in the module doesn't really matter
        