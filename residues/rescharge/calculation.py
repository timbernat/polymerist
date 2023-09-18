'''For determining library charges from a molecule which with assigned partial charges'''

import logging
LOGGER = logging.getLogger(__name__)

from typing import Optional
from abc import ABC, abstractmethod
from dataclasses import dataclass

from collections import defaultdict
from pathlib import Path
import numpy as np

from rdkit import Chem
from rdkit.Chem import Mol as RDMol

from openff.toolkit.topology.molecule import Molecule
from openmm.unit import elementary_charge

from .chgtypes import ChargesByResidue, ChargeMap
from ...genutils.maths.statistics import Accumulator
from ...monomers.repr import MonomerGroup


# INTERFACE FOR DISTRIBUTING EXCESS RESIDUE CHARGES
@dataclass
class ChargeRedistributionStrategy(ABC):
    '''Interface for defining how any excess charge should be distributed within residues to ensure a given overall net charge'''
    desired_net_charge : float = 0.0 # by default, make neutral

    @abstractmethod # !NOTE! : this must be implemented in child classes
    def _determine_charge_offsets(self, base_charges : ChargeMap, fragment : RDMol, net_charge_diff : float) -> ChargeMap:
        '''Provided a set of base charges, a structural molecule fragment, and a desired net charge,
        determine what charges offsets need to be applied where in order to achieve the desired net charge'''
        raise NotImplemented

    def redistributed_charges(self, base_charges : ChargeMap, fragment : RDMol) -> ChargeMap:
        '''Take a map of base charges and a structural fragment for a residue and a desired net charge (typically neutral, i.e. 0)
        and return a new charge map with the excess/deficit charge distributed in such a way as to make the residue have the desired net charge'''
        net_charge_diff = self.desired_net_charge - sum(chg for chg in base_charges.values())
        charge_offsets = self._determine_charge_offsets(base_charges, fragment, net_charge_diff)

        new_charges = {
            sub_id : charge + charge_offsets[sub_id]
                for sub_id, charge in base_charges.items()
        }

        # assert(sum(chg for chg in new_charges.values()) == desired_net_charge) # double check charges were correctly redistributed - TODO : find more reliable way to check this than floating-point comparison
        return new_charges
            
class UniformDistributionStrategy(ChargeRedistributionStrategy):
    '''Simplest possible strategy, distribute any excess charge in a residue according to a uniform distribution (spread evenly)'''
    def _determine_charge_offsets(base_charges : ChargesByResidue, fragment : RDMol, net_charge_diff : float) -> ChargesByResidue:
        return {sub_id : net_charge_diff / len(base_charges) for sub_id in base_charges}


# functions for determining library charges
def find_repr_residues(mol : Molecule) -> dict[str, int]:
    '''Determine names and smallest residue numbers of all unique residues in charged molecule
    Used as representatives for generating labelled SMARTS strings '''
    rep_res_nums = defaultdict(set) # numbers of representative groups for each unique residue, used to build SMARTS strings
    for atom in mol.atoms: 
        rep_res_nums[atom.metadata['residue_name']].add(atom.metadata['residue_number']) # collect unique residue numbers

    for res_name, ids in rep_res_nums.items():
        rep_res_nums[res_name] = min(ids) # choose group with smallest id of each residue to denote representative group

    return rep_res_nums

def get_averaged_charges(cmol : Molecule, monomer_group : MonomerGroup, cds : Optional[ChargeRedistributionStrategy]=UniformDistributionStrategy(desired_net_charge=0.0)) -> ChargesByResidue:
    '''Takes a charged molecule and a dict of monomer SMIRKS strings and averages charges for each repeating residue. 
    Returns a list of ChargedResidue objects, each of which holds:
        - A dict of the averaged charges by atom 
        - The name of the residue associated with the charges
        - A SMARTS string of the residue's structure
        - An nx.Graph representing the structure of the residue'''
    # rdmol = cmol.to_rdkit() # create rdkit representation of Molecule to allow for SMARTS generation
    rep_res_nums = find_repr_residues(cmol) # determine ids of representatives of each unique residue

    atom_id_mapping   = defaultdict(lambda : defaultdict(int))
    res_charge_accums = defaultdict(lambda : defaultdict(Accumulator))
    for atom in cmol.atoms: # accumulate counts and charge values across matching subsftructures
        res_name, res_num     = atom.metadata['residue_name'   ], atom.metadata['residue_number']
        substruct_id, atom_id = atom.metadata['substructure_id'], atom.metadata['pdb_atom_id'   ]

        if res_num == rep_res_nums[res_name]: # if atom is member of representative group for any residue...
            # rdmol.GetAtomWithIdx(atom_id).SetAtomMapNum(atom_id)  # ...and set atom number for labelling in SMARTS string
            atom_id_mapping[res_name][atom_id] = (substruct_id, atom.symbol) # ...collect pdb id...

        curr_accum = res_charge_accums[res_name][substruct_id] # accumulate charge info for averaging
        curr_accum.sum += atom.partial_charge.magnitude # eschew units (easier to handle, added back when writing to XML)
        curr_accum.count += 1

    chgs_by_res = ChargesByResidue
    for res_name, charge_accums in res_charge_accums.items():
        # rdSMARTS = rdmolfiles.MolFragmentToSmarts(rdmol, atomsToUse=atom_id_mapping[res_name].keys()) # determine SMARTS for the current residue's representative group
        # mol_frag = rdmolfiles.MolFromSmarts(rdSMARTS) # create fragment from rdkit SMARTS to avoid wild atoms (using rdkit over nx.subgraph for more detailed atomwise info)
        
        SMARTS = monomer_group.monomers[res_name] # extract SMARTS string from monomer data

        charge_map = {substruct_id : accum.average for substruct_id, accum in charge_accums.items()} 
        if cds is not None:
            mol_frag = Chem.MolFromSmarts(SMARTS) # TODO : make this take specific fragments from the original molecule in an id-preserving way
            charge_map = cds.redistributed_charges(charge_map, mol_frag)

        chgs_by_res.charges[res_name] = charge_map

    return chgs_by_res