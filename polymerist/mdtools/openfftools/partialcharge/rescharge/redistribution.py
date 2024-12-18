'''Strategies for redistribution excess partial charge among residues'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from dataclasses import dataclass
from abc import ABC, abstractmethod
from rdkit.Chem import Mol

from .rctypes import ChargeMap, ChargesByResidue


# ABSTRACT INTERFACE FOR DISTRIBUTING EXCESS RESIDUE CHARGES
@dataclass
class ChargeRedistributionStrategy(ABC):
    '''Interface for defining how any excess charge should be distributed within residues to ensure a given overall net charge'''
    desired_net_charge : float = 0.0 # by default, make neutral

    @abstractmethod # !NOTE! : this must be implemented in child classes
    def _determine_charge_offsets(self, base_charges : ChargeMap, fragment : Mol, net_charge_diff : float) -> ChargeMap:
        '''Provided a set of base charges, a structural molecule fragment, and a desired net charge,
        determine what charges offsets need to be applied where in order to achieve the desired net charge'''
        raise NotImplementedError

    def redistributed_charges(self, base_charges : ChargeMap, fragment : Mol) -> ChargeMap:
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

## CONCRETE IMPLEMENTATIONS OF REDISTRIBUTION STRATEGIES
class UniformDistributionStrategy(ChargeRedistributionStrategy):
    '''Simplest possible strategy, distribute any excess charge in a residue according to a uniform distribution (spread evenly)'''
    def _determine_charge_offsets(self, base_charges : ChargesByResidue, fragment : Mol, net_charge_diff : float) -> ChargesByResidue:
        return {sub_id : net_charge_diff / len(base_charges) for sub_id in base_charges}
