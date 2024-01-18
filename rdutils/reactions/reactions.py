'''Classes for representing information about reaction mechanisms and tracing bonds and atoms along a reaction'''

from typing import Iterable, Optional, Union
from dataclasses import dataclass, field

from io import StringIO
from pathlib import Path
from functools import cached_property

from rdkit import Chem
from rdkit.Chem import rdChemReactions

from ..rdtypes import RDMol
from ..labeling.molwise import ordered_map_nums
from ..labeling.bondwise import get_bonded_pairs_by_map_nums
from ..amalgamation import bonding

from ...genutils.fileutils.pathutils import aspath, asstrpath


# REACTION INFORMATICS CLASSES
@dataclass
class RxnProductInfo:
    '''For storing atom map numbers associated with product atoms and bonds participating in a reaction'''
    prod_num : int
    reactive_atom_map_nums : list[int] = field(default_factory=list)

    new_bond_ids_to_map_nums : dict[int, tuple[int, int]] = field(default_factory=dict)
    mod_bond_ids_to_map_nums : dict[int, tuple[int, int]] = field(default_factory=dict)
    
class AnnotatedReaction(rdChemReactions.ChemicalReaction):
    '''
    RDKit ChemicalReaction with additional useful information about product atom and bond mappings

    Initialization must be done either via AnnotatedReaction.from_smarts or AnnotatedReaction.from_rdmols,
    asdirect override of pickling in __init__ method not yet implemented
    '''

    # LOADING METHODS
    @classmethod
    def from_smarts(cls, rxn_smarts : str) -> 'AnnotatedReaction':
        '''For instantiating reactions from SMARTS strings'''
        return cls(rdChemReactions.ReactionFromSmarts(rxn_smarts)) # pass to init method

    @classmethod
    def from_rxnfile(cls, rxnfile_path : Union[str, Path]) -> 'AnnotatedReaction':
        '''For instantiating reactions directly from MDL .rxn files'''
        return cls(rdChemReactions.ReactionFromRxnFile(asstrpath(rxnfile_path)))

    @classmethod
    def from_rdmols(cls, reactant_templates : Iterable[RDMol], product_templates : Iterable[RDMol]) -> 'AnnotatedReaction':
        '''For instantiating reactions directly from molecules instead of SMARTS strings'''
        # 1) label atoms as belonging to reactant or product via RDKit 'magic' internal property (1 = reactant, 2 = product, 3 = agent)
        for reactant in reactant_templates: # TODO : implement non-in-place assignment of these properties
            for atom in reactant.GetAtoms():
                atom.SetIntProp('molRxnRole', 1) 

        for product in product_templates: # TODO : implement non-in-place assignment of these properties
            for atom in product.GetAtoms():
                atom.SetIntProp('molRxnRole', 2) 

        # 2) generate single combined molecule to fit RDKit reaction spec
        rxn_mol = bonding.combined_rdmol(*reactant_templates, *product_templates, assign_map_nums=False, editable=False) # kwargs are explicitly needed here

        return cls(rdChemReactions.ReactionFromMolecule(rxn_mol))
    
    # I/O METHODS
    def to_rxnfile(self, rxnfile_path : Union[str, Path], rxnname : Optional[str]=None, _rxnname_line : int=1) -> None:
        '''Save reaction to an MDL .RXN file. Replaces ports with R-groups to enable proper loading'''
        rxn_block = rdChemReactions.ReactionToRxnBlock(self)
        rxn_block = rxn_block.replace('*', 'R')

        with aspath(rxnfile_path).open('w') as rxnfile: # TODO : replace aspath with allow_str_paths decorator
            for i, line in enumerate(StringIO(rxn_block)):
                if (rxnname is not None) and (i == _rxnname_line): # name inserted into second line per CTFile spec (https://discover.3ds.com/sites/default/files/2020-08/biovia_ctfileformats_2020.pdf) 
                    rxnfile.write(f'\t\t\t{rxnname}\n')
                else:
                    rxnfile.write(line)

    # RXN INFO METHODS
    @cached_property
    def reacting_atom_map_nums(self) -> list[int]:
        '''List of the map numbers of all reactant atoms which participate in the reaction'''
        return [
            self.GetReactantTemplate(reactant_id).GetAtomWithIdx(atom_id).GetAtomMapNum()
                for reactant_id, reacting_atom_ids in enumerate(self.GetReactingAtoms())
                    for atom_id in reacting_atom_ids
        ]
    
    @cached_property
    def map_nums_to_reactant_nums(self) -> dict[int, int]:
        '''Back-map yielding the index of the source reactant for the atom of each map number'''
        return {
            atom.GetAtomMapNum() : i
                for i, react_template in enumerate(self.GetReactants())
                    for atom in react_template.GetAtoms()
        }
    
    @cached_property
    def product_info_maps(self) -> dict[int, RxnProductInfo]:
        '''Map from product index to information about reactive atoms and bonds in that product'''
        # map reacting atoms and bonds for each product
        prod_info_map = {}
        for i, product_template in enumerate(self.GetProducts()):
            prod_info = RxnProductInfo(i)
            prod_info.reactive_atom_map_nums = [
                map_num
                    for map_num in self.reacting_atom_map_nums
                        if map_num in ordered_map_nums(product_template)
            ]

            for bond_id, atom_id_pair in get_bonded_pairs_by_map_nums(product_template, *prod_info.reactive_atom_map_nums).items(): # consider each pair of reactive atoms
                map_num_1, map_num_2 = map_num_pair = tuple(
                    product_template.GetAtomWithIdx(atom_id).GetAtomMapNum()
                        for atom_id in atom_id_pair
                )
                
                if self.map_nums_to_reactant_nums[map_num_1] == self.map_nums_to_reactant_nums[map_num_2]: # if reactant IDs across bond match, the bond must have been modified (i.e. both from single reactant...)
                    prod_info.mod_bond_ids_to_map_nums[bond_id] = map_num_pair
                else: # otherwise, bond must be newly formed (spans between previously disjoint monomers) 
                    prod_info.new_bond_ids_to_map_nums[bond_id] = map_num_pair
            prod_info_map[i] = prod_info
        
        return prod_info_map
    
