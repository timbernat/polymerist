'''Classes for representing information about reaction mechanisms and tracing bonds and atoms along a reaction'''

from typing import ClassVar, Iterable, Optional, Sequence, Union
from dataclasses import dataclass, field

import re
from io import StringIO
from pathlib import Path
from functools import cached_property

from rdkit.Chem import rdChemReactions

from .reactexc import BadNumberReactants
from ..rdtypes import RDMol
from ..bonding._bonding import combined_rdmol
from ..labeling.molwise import ordered_map_nums
from ..labeling.bondwise import get_bonded_pairs_by_map_nums

from ...genutils.decorators.functional import allow_string_paths, allow_pathlib_paths
from ...genutils.sequences import bin_ids_forming_sequence
from ...smileslib.substructures import matching_labels_from_substruct_dict


# REACTION INFORMATICS CLASSES
@dataclass
class RxnProductInfo:
    '''For storing atom map numbers associated with product atoms and bonds participating in a reaction'''
    prod_num : int
    reactive_atom_map_nums : list[int] = field(default_factory=list)

    new_bond_ids_to_map_nums : dict[int, tuple[int, int]] = field(default_factory=dict)
    mod_bond_ids_to_map_nums : dict[int, tuple[int, int]] = field(default_factory=dict)
    
RXNNAME_RE = re.compile(r'^\t*(?P<rxnname>.*?)\n$')
class AnnotatedReaction(rdChemReactions.ChemicalReaction):
    '''
    RDKit ChemicalReaction subclass with additional useful information about product atom and bond mappings and reaction naming
    Initialization must be done either via AnnotatedReaction.from_smarts, AnnotatedReaction.from_rdmols, or AnnotatedReaction.from_rxnfile
    '''
    # line number in .rxn file where (optional) name of reaction should be located (per the CTFile spec https://discover.3ds.com/sites/default/files/2020-08/biovia_ctfileformats_2020.pdf) 
    _RXNNAME_LINE_NO : ClassVar[int] = 1
    
    # LOADING METHODS
    @classmethod
    def from_smarts(cls, rxn_smarts : str) -> 'AnnotatedReaction':
        '''For instantiating reactions from SMARTS strings'''
        return cls(rdChemReactions.ReactionFromSmarts(rxn_smarts)) # pass to init method

    @classmethod
    def from_rdmols(cls, reactant_templates : Iterable[RDMol], product_templates : Iterable[RDMol], agent_templates : Optional[Iterable[RDMol]]=None) -> 'AnnotatedReaction':
        '''For instantiating reactions directly from molecules instead of SMARTS strings'''
        # label atoms as belonging to reactant or product via RDKit 'magic' internal property (1 = reactant, 2 = product, 3 = agent)
        if agent_templates is None:
            agent_templates = []

        template_role_map = { # RDKit "magic" aliases which define what role each atom plays in a reaction
            1 : reactant_templates,
            2 : product_templates,
            3 : agent_templates,
        }

        for mol_rxn_role, templates in template_role_map.items():
            for template in templates: # TODO : implement non-in-place assignment of these properties
                for atom in template.GetAtoms():
                    atom.SetIntProp('molRxnRole', mol_rxn_role) 
        rxn_mol = combined_rdmol(*reactant_templates, *product_templates, *agent_templates, assign_map_nums=False, editable=False) # kwargs are explicitly needed here

        return cls(rdChemReactions.ReactionFromMolecule(rxn_mol))

    # I/O METHODS
    @classmethod
    @allow_string_paths
    def rxnname_from_rxnfile(cls, rxnfile_path : Union[str, Path]) -> str:
        '''Extract the reaction name from a (properly-formatted) .RXN file'''
        with rxnfile_path.open('r') as file:
            for i, line in enumerate(file):
                if i == cls._RXNNAME_LINE_NO:
                    return re.match(RXNNAME_RE, line).group('rxnname')
                
    @property
    def rxnname(self) -> str:
        '''A string handle associated with this reaction'''
        return getattr(self, '_rxnname', '') # default to the empty string if unset
    
    @rxnname.setter
    def rxnname(self, new_rxnname : str) -> None:
        '''Set a new name for the current reaction'''
        self._rxnname = new_rxnname

    @classmethod
    @allow_pathlib_paths
    def from_rxnfile(cls, rxnfile_path : Union[str, Path]) -> 'AnnotatedReaction':
        '''For instantiating reactions directly from MDL .rxn files'''
        rxn = cls(rdChemReactions.ReactionFromRxnFile(rxnfile_path))
        rxn.rxnname = cls.rxnname_from_rxnfile(rxnfile_path)

        return rxn
    
    @allow_string_paths
    def to_rxnfile(self, rxnfile_path : Union[str, Path], wilds_to_R_groups : bool=True) -> None:
        '''Save reaction to an MDL .RXN file. Replaces ports with R-groups to enable proper loading'''
        rxn_block = rdChemReactions.ReactionToRxnBlock(self)
        if wilds_to_R_groups:
            rxn_block = rxn_block.replace('*', 'R')

        with rxnfile_path.open('w') as rxnfile:
            for i, line in enumerate(StringIO(rxn_block)):
                if (i == self._RXNNAME_LINE_NO):
                    rxnfile.write(f'\t\t\t{self.rxnname}\n')
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
    
    def valid_reactant_ordering(self, reactants : Sequence[RDMol]) -> Optional[list[RDMol]]:
        '''Given an RDKit chemical reaction mechanism and a sequence of reactant Mols, will determine if there is
        an ordering of the reactants which is compatible with the reactant templates defined in the reaction

        Returns the first found ordering, or NoneType if no such ordering exists'''
        reactant_templates_by_index = {i : reac_templ for i, reac_templ in enumerate(self.GetReactants())}
        num_reactants_in_mechanism = self.GetNumReactantTemplates() # = len(reactant_templates_by_index)
        num_reactants_provided = len(reactants)

        # preliminary (quick) check; are there the right number of Mols
        if num_reactants_provided != num_reactants_in_mechanism:
            raise BadNumberReactants(f'{self.__class__.__name__} expected {num_reactants_in_mechanism} reactants, but {num_reactants_provided} were provided')

        # if number of fragments is correct, perform more complex subset choice evaluation
        possible_fragment_orderings = bin_ids_forming_sequence( # generates all possible orders of the fragments which match the expected reactant templates for the reaction
            sequence=reactant_templates_by_index.keys(),
            choice_bins = ( # generator (rather than list) comprehension will suffice here, since there is no need to reuse these bins
                matching_labels_from_substruct_dict(reactant, reactant_templates_by_index)
                    for reactant in reactants
            ),
            unique_bins=True
        )
        try:
            return [reactants[i] for i in next(possible_fragment_orderings)] # get first valid ordering
        except StopIteration:
            return None # slightly verbose, but prefer to be explicit here