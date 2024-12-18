'''Classes for representing information about reaction mechanisms and tracing bonds and atoms along a reaction'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import ClassVar, Iterable, Optional, Sequence, Union
from dataclasses import dataclass, field

import re
from io import StringIO
from pathlib import Path
from functools import cached_property

from rdkit.Chem import rdChemReactions, Mol

from .reactexc import BadNumberReactants
from ..bonding._bonding import combined_rdmol
from ..labeling.molwise import ordered_map_nums
from ..labeling.bondwise import get_bonded_pairs_by_map_nums

from ...genutils.decorators.functional import allow_string_paths, allow_pathlib_paths
from ...genutils.sequences.discernment import DISCERNMENTSolver
from ...smileslib.substructures import num_substruct_queries_distinct


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
    
    # LOADING/EXPORT METHODS
    @classmethod
    def from_smarts(cls, rxn_smarts : str) -> 'AnnotatedReaction':
        '''Iinstantiate reaction from mapped SMARTS string'''
        return cls(rdChemReactions.ReactionFromSmarts(rxn_smarts))
    
    # NOTE : cannot analogous implement "from_smiles" classmethod, as rdChemreactions does not support initialization from SMILES (only SMARTS)

    def to_smarts(self) -> str:
        '''Export reaction as mapped SMARTS string'''
        return rdChemReactions.ReactionToSmarts(self) # TODO : implement * -> R replacement here (rather than in rxn file I/O)

    @property
    def smarts(self) -> str:
        '''Mapped SMARTS string representation of reaction'''
        return self.to_smarts()

    def to_smiles(self) -> str:
        '''Export reaction as mapped SMILES string'''
        return rdChemReactions.ReactionToSmiles(self) # TODO : implement * -> R replacement here (rather than in rxn file I/O)

    @property
    def smiles(self) -> str:
        '''Mapped SMILES string representation of reaction'''
        return self.to_smiles()

    @classmethod
    def from_rdmols(cls, reactant_templates : Iterable[Mol], product_templates : Iterable[Mol], agent_templates : Optional[Iterable[Mol]]=None) -> 'AnnotatedReaction':
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
    
    def valid_reactant_ordering(self, reactants : Sequence[Mol], as_mols : bool=True) -> Optional[list[Mol]]:
        '''
        Given an RDKit chemical reaction mechanism and a sequence of reactant Mols, will determine if there is
        an ordering of the reactants which is compatible with the reactant templates defined in the reaction

        Returns the first found ordering, or NoneType if no such ordering exists
        Ordering returned as list of Chem.Mol objects if as_mols == True, or as list of ints otherwise
        '''
        # 0) Preliminary quick check on number of reactants; can discount a bad reactant collection prior to more expensive check
        num_reactants_provided = len(reactants)
        num_reactant_templates_in_mechanism = self.GetNumReactantTemplates()

        if num_reactants_provided != num_reactant_templates_in_mechanism:
            raise BadNumberReactants(f'{self.__class__.__name__} expected {num_reactant_templates_in_mechanism} reactants, but {num_reactants_provided} were provided')

        # if number of fragments is correct, perform more complex of whether a molecule-unique substructure selection exists
        ## 1) Pre-compilation of template substructure indices
        reactant_template_indices : list[int] = [templ_idx for templ_idx in range(num_reactant_templates_in_mechanism)]
        reactant_templates_by_index : dict[int, Mol] = {
            templ_idx : self.GetReactantTemplate(templ_idx)
                for templ_idx in reactant_template_indices
        }

        template_substruct_ids_by_reactant_ids : dict[int, list[int]] = { 
            react_idx : [ # associate with each reactant a substructure "word"...
                templ_idx # ... consisting of the indices of reactant template indices present as substructure(s) in the reactant...
                    for templ_idx, react_templ in reactant_templates_by_index.items()
                        for _ in range(num_substruct_queries_distinct(reactant, react_templ)) # ...which appear as many times as the substructure is uniquely present
            ]
                for react_idx, reactant in enumerate(reactants)
        } # TOSELF: implemented this w/ indices (rather than directly using Mols, which are hashable) because the same mol might return distinct hashes (cannot evaluate self-equality via __eq__ directly)

        ## 2) Choice-finding to see if a selection of reactants matching the templates exists
        reactant_ordering_planner = DISCERNMENTSolver(template_substruct_ids_by_reactant_ids)

        if reactant_ordering_planner.solution_exists(reactant_template_indices, unique_bins=True):
            reactant_ordering : tuple[int] = next(reactant_ordering_planner.enumerate_choices(
                word=reactant_template_indices,
                unique_bins=True) # need to have unique bins so a single reactant is not taken more than once
            )
            if as_mols:
                return [reactants[i] for i in reactant_ordering] # get first valid ordering
            return reactant_ordering
        else: # NOTE: this else clause is not strictly necessary (would otherwise return None anyway), but prefer to have it for explicitness
            return None 