'''Classes for representing information about reaction mechanisms and tracing bonds and atoms along a reaction'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import ClassVar, Generator, Iterable, Optional, Sequence, Union
from dataclasses import dataclass, field
from enum import StrEnum, auto

import re
from io import StringIO
from pathlib import Path

from functools import cached_property
from collections import defaultdict, Counter

from rdkit.Chem import rdChemReactions, Mol, BondType

from ..chemselect import get_mapped_atoms
from ..bonding._bonding import combined_rdmol
from ..labeling.molwise import ordered_map_nums
from ..labeling.bondwise import get_bonded_pairs_by_map_nums

from ...smileslib.sanitization import canonical_SMILES_from_mol
from ...smileslib.substructures import num_substruct_queries_distinct
from ...genutils.decorators.functional import allow_string_paths, allow_pathlib_paths
from ...genutils.sequences.discernment import DISCERNMENTSolver, SymbolInventory


# REACTION INFORMATICS CLASSES
class BondChange(StrEnum):
    '''For indicating how a bond which changed in a reaction was altered'''
    ADDED = auto()
    DELETED = auto()
    MODIFIED = auto() # specifically, when bond order is modified but the bond persists
    UNCHANGED = auto()

@dataclass(frozen=True)
class AtomTraceInfo:
    '''For encapsulating information about the origin and destination of a mapped atom, traced through a reaction'''
    map_num : int
    reactant_idx      : int # index of the reactant template within a reaction in which the atom occurs
    reactant_atom_idx : int # index of the target atom WITHIN the above reactant template
    product_idx       : int # index of the product template within a reaction in which the atom occurs
    product_atom_idx  : int # index of the target atom WITHIN the above product template
   
@dataclass(frozen=True)
class BondTraceInfo:
    '''For encapsulating information about bonds which are between mapped atoms and which change during a reaction'''
    map_nums : tuple[int, int] # map numbers of the pair of atoms the bond connects
    # NOTE: reactant index doesn't make much sense, since the atoms the bond spans might have comes from two distinct reactant templates
    product_idx      : int # index of the reactant template within a reaction in which the modified bond occurs
    product_bond_idx : int # index of the target bond WITHIN the above product template
    change_type  : Union[str, BondChange]
    bond_order   : Union[float, BondType] # bond order in the product (i.e. AFTER the change)
    

@dataclass
class RxnProductInfo:
    '''For storing atom map numbers associated with product atoms and bonds participating in a reaction'''
    prod_num : int
    reactive_atom_map_nums : list[int] = field(default_factory=list)

    new_bond_ids_to_map_nums : dict[int, tuple[int, int]] = field(default_factory=dict)
    mod_bond_ids_to_map_nums : dict[int, tuple[int, int]] = field(default_factory=dict)
    
    
# REACTION CLASS
RXNNAME_RE = re.compile(r'^\t*(?P<rxnname>.*?)\n$')
class AnnotatedReaction(rdChemReactions.ChemicalReaction):
    '''
    RDKit ChemicalReaction subclass with additional useful information about product atom and bond mappings and reaction naming
    Initialization must be done either via AnnotatedReaction.from_smarts, AnnotatedReaction.from_rdmols, or AnnotatedReaction.from_rxnfile
    '''
    # line number in .rxn file where (optional) name of reaction should be located (per the CTFile spec https://discover.3ds.com/sites/default/files/2020-08/biovia_ctfileformats_2020.pdf) 
    _RXNNAME_LINE_NO : ClassVar[int] = 1
    _atom_ridx_prop_name   : ClassVar[str] = 'reactant_idx' # name of the property to assign reactant indices to; set for entire class
    _bond_change_prop_name : ClassVar[str] = 'bond_changed' # name of property to set on bonds to indicated they have changed in a reaction

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)

    # LOADING/EXPORT METHODS
    @classmethod
    def from_smarts(cls, rxn_smarts : str) -> 'AnnotatedReaction':
        '''Instantiate reaction from mapped SMARTS string'''
        return cls(rdChemReactions.ReactionFromSmarts(rxn_smarts.replace('#0', '*'))) # clean up any wild-atom conversion artifacts from porting a SMARTS through SMILES
    
    # NOTE : cannot analogous implement "from_smiles" classmethod, as rdChemreactions does not support initialization from SMILES (only SMARTS)

    def to_smarts(self) -> str:
        '''Export reaction as mapped SMARTS string''' # TODO : implement * -> R replacement here (rather than in rxn file I/O)
        return rdChemReactions.ReactionToSmarts(self).replace('#0', '*') # clean up any wild-atom conversion artifacts from porting a SMARTS through SMILES 

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

    # INFO FOR TRACING ATOMS THROUGH REACTIONS
    @cached_property
    def reactant_idxs_by_map_num(self) -> dict[int, tuple[int, int]]:
        '''
        Associates map numbers of all mapped atoms to the pair of reactant index and
        the atom index with that reactant, respectively, that the mapped atom apears in 
        
        Keyed by map number, with values of (reactant index, atom index) tuples
        '''
        return {
            atom.GetAtomMapNum() : (reactant_template_idx, atom.GetIdx())
                for reactant_template_idx, reactant_template in enumerate(self.GetReactants())
                    for atom in get_mapped_atoms(reactant_template, as_indices=False)
        }

    @cached_property
    def product_idxs_by_map_num(self) -> dict[int, tuple[int, int]]:
        '''
        Associates map numbers of all mapped atoms to the pair of product index and
        the atom index with that product, respectively, that the mapped atom apears in 
        
        Keyed by map number, with values of (product index, atom index) tuples
        '''
        return {
            atom.GetAtomMapNum() : (product_template_idx, atom.GetIdx())
                for product_template_idx, product_template in enumerate(self.GetProducts())
                    for atom in get_mapped_atoms(product_template, as_indices=False)
        }
    
    @cached_property
    def mapped_atom_info(self) -> dict[int, AtomTraceInfo]:
        '''Compile reactant origin and product destination of all mapped atoms'''
        mapped_atom_infos = {}
        for map_number, (reactant_idx, reactant_atom_idx) in self.reactant_idxs_by_map_num.items():
            product_idx, product_atom_idx = self.product_idxs_by_map_num[map_number]
            mapped_atom_infos[map_number] = AtomTraceInfo(
                map_num=map_number, # DEVNOTE: am aware this information seems redundant, but this duplication allows one to vary the level of coupling to the returned dict
                reactant_idx=reactant_idx,
                reactant_atom_idx=reactant_atom_idx,
                product_idx=product_idx,
                product_atom_idx=product_atom_idx,
            )
        return mapped_atom_infos
           
    @cached_property
    def reactive_atom_info(self) -> dict[int, AtomTraceInfo]:
        '''Compile reactant origin and product destination of all mapped atoms which are changed by the reaction'''
        reactive_atom_infos = {}
        for reactant_idx, reactant_atom_idxs in enumerate(self.GetReactingAtoms(mappedAtomsOnly=True)):
            reactant_template = self.GetReactantTemplate(reactant_idx)
            for reactant_atom_idx in reactant_atom_idxs:
                map_num = reactant_template.GetAtomWithIdx(reactant_atom_idx).GetAtomMapNum()
                reactive_atom_infos[map_num] = self.mapped_atom_info[map_num]
        
        return reactive_atom_infos
        
    @cached_property
    def reacting_atom_map_nums(self) -> list[int]:
        '''List of the map numbers of all reactant atoms which participate in the reaction'''
        return [
            self.GetReactantTemplate(reactant_id).GetAtomWithIdx(atom_id).GetAtomMapNum()
                for reactant_id, reacting_atom_ids in enumerate(self.GetReactingAtoms())
                    for atom_id in reacting_atom_ids
        ]
        
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
                
                if self.reactant_idxs_by_map_num[map_num_1][0] == self.reactant_idxs_by_map_num[map_num_2][0]: # if reactant IDs across bond match, the bond must have been modified (i.e. both from single reactant...)
                    prod_info.mod_bond_ids_to_map_nums[bond_id] = map_num_pair
                else: # otherwise, bond must be newly formed (spans between previously disjoint monomers) 
                    prod_info.new_bond_ids_to_map_nums[bond_id] = map_num_pair
            prod_info_map[i] = prod_info
        
        return prod_info_map
    
    # REACTANT PATHFINDING
    def compile_functional_group_inventory(self, reactants : Iterable[Mol], label_reactants_with_smiles : bool=False) -> SymbolInventory:
        '''
        Construct an inventory of numbers of functional groups (reactant templates) found in a sequence of reactant Mols,
        which can be evaluated as a DISCERNMENT-type problem to determine valid reactant ordering(s), if some exist
        
        Isomorphisms to the DISCERNMENT problem in this instance are as follows:
            "symbols"     <-> indices of functional groups, as-defined by reactant templates
            "target word" <-> the sequence of indices [n-1] = [0, 1, ..., n-1] where "n" is the number of reactant templates
            "magazine"    <-> an ordered collection of reactant molecules, which may contain any number of functional groups apiece each
            "word labels" <-> either the index of a reactant in the order it's passed or its canonical SMILES representation (depends on value of "label_reactants_with_smiles")
        '''
        fn_group_sym_inv = defaultdict(Counter) # TODO: add more "native" mechanism to instantiate SymbolInventory directly from counts (rather than sequences of bins)
        for reactant_idx, reactant in enumerate(reactants): # NOTE: placed first in case reactants are exhaustible (e.g. generator-like)
            for template_idx, reactant_template in enumerate(self.GetReactants()): 
                fn_group_sym_inv[template_idx][
                    canonical_SMILES_from_mol(reactant) if label_reactants_with_smiles else reactant_idx
                ] = num_substruct_queries_distinct(reactant, reactant_template) # counts are the number of times a particular reactant template occurs in a monomer
        
        return SymbolInventory(fn_group_sym_inv)
        
    def enumerate_valid_reactant_orderings(
            self,
            reactants : Sequence[Mol],
            as_mols : bool=True,
            allow_resampling : bool=False
        ) -> Generator[Union[None, tuple[int], tuple[Mol]], None, None]:
        '''
        Enumerates all orderings of reactants compatible with the reactant templates defined in this reaction

        Yields:
        * a single NoneType if no such ordering exists
        * tuples of Chem.Mol objects if as_mols=True
        * tuples of indices of reactants in the passed sequence if as_mols=False
        
        If allow_resampling=False, each reactant will only be allowed to contribute exactly 1 of its functional groups 1 to any solution 
        '''
        reactant_orderings_found : bool = False
        reactant_ordering_planner = DISCERNMENTSolver(self.compile_functional_group_inventory(reactants, label_reactants_with_smiles=False)) 
        for reactant_ordering in reactant_ordering_planner.enumerate_choices(
            word=[templ_idx for templ_idx in range(self.GetNumReactantTemplates())], # use indices of reactants to allow direct lookup of reactants from solution
            unique_bins=not allow_resampling, # if one is OK with the same reactant being used for multiple functional groups, ignore its functional group multiplicities
        ):                                          
            yield tuple(reactants[i] for i in reactant_ordering) if as_mols else tuple(reactant_ordering)
            reactant_orderings_found = True # flag that at least one ordering is possible
        
        if not reactant_orderings_found:  # if no solution exists, explicitly yield a single NoneType sentinel value (simplifies check for no solutions existing)
            yield None 
    
    def valid_reactant_ordering(
            self,
            reactants : Sequence[Mol],
            as_mols : bool=True,
            allow_resampling : bool=False
        ) -> Union[None, tuple[int], tuple[Mol]]:
        '''
        Get first ordering of reactants compatible with the reactant templates defined in this reaction

        Yields:
        * a single NoneType if no such ordering exists
        * a tuple of Chem.Mol objects if as_mols=True
        * a tuple of indices of reactants in the passed sequence if as_mols=False
        
        If allow_resampling=False, each reactant will only be allowed to contribute exactly 1 of its functional groups 1 to any solution 
        '''
        return next(self.enumerate_valid_reactant_orderings(
            reactants,
            as_mols=as_mols,
            allow_resampling=allow_resampling
        )) # DEVNOTE: enumeration guarantees a NoneType is returned when no solution exists (no edge-case handling needed!)
        
    def has_reactable_subset(self, reactants : Sequence[Mol], allow_resampling : bool=False) -> bool:
        '''
        Determine if a sequence of reactants Mols contains any subset of Mols which are compatible with the reactant templates defined by this reaction
        If allow_resampling=False, each reactant will only be allowed to contribute exactly 1 of its functional groups 1 to any solution 
        '''
        return self.valid_reactant_ordering(reactants=reactants, as_mols=False, allow_resampling=allow_resampling) is not None
