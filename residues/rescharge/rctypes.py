'''Custom types used in determining residue charges'''

from typing import TypeAlias
from dataclasses import dataclass, field

from rdkit.Chem import Mol as RDMol
from ...genutils.fileutils.jsonio.jsonify import make_jsonifiable


ChargeMap : TypeAlias = dict[int, float] # maps substructure IDs to partial charge values

@dataclass
class ChargedResidue:
    '''Dataclass for more conveniently storing averaged charges for a residue group'''
    charges : ChargeMap
    residue_name : str
    SMARTS : str
    mol_fragment : RDMol

@make_jsonifiable
@dataclass
class ChargesByResidue:
    '''Class for storing substructure charge maps by residue'''
    charges : dict[str, ChargeMap] = field(default_factory=dict)
