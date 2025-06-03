'''Classes for representing information about reaction mechanisms and tracing bonds and atoms along a reaction'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import (
    Callable,
    ClassVar,
    Generator,
    Iterable,
    Mapping,
    Optional,
    Sequence,
    TypeVar,
    Union,
)
L = TypeVar('L') # generic type for molecule labels

import re
from io import StringIO
from pathlib import Path

from random import shuffle
from itertools import chain
from functools import cached_property
from collections import defaultdict, Counter

from rdkit.Chem.rdchem import Mol, Atom, Bond, BondType
from rdkit.Chem.rdChemReactions import (
    ChemicalReaction,
    ReactionFromSmarts,
    ReactionToSmarts,
    ReactionToSmiles,
    ReactionFromRxnBlock,
    ReactionFromRxnFile,
    ReactionToRxnBlock,
    ReactionFromMolecule,
)
from .reactexc import BadNumberReactants, ReactantTemplateMismatch
from .reactinfo import AtomTraceInfo, BondTraceInfo, BondChange, REACTANT_INDEX_PROPNAME, BOND_CHANGE_PROPNAME

from ..rdprops import copy_rdobj_props
from ..bonding import combined_rdmol
from ..selection import mapped_atoms, mapped_neighbors
from ..sanitization import sanitizable_mol_outputs
from ..substructures import num_substruct_queries_distinct

from ...genutils.fileutils.pathutils import allow_string_paths, allow_pathlib_paths
from ...genutils.sequences.discernment import DISCERNMENTSolver, SymbolInventory


# helper functions
def map_numbers_to_neighbor_bonds(mol : Mol, atom_idx : int) -> dict[int, BondType]:
    '''
    Given an atom, get a mapping from the atom map numbers of explicitly-mapped
    neighbor atoms to the bond types of the bonds connecting the atoms
    '''
    return {
        neighbor_atom.GetAtomMapNum() : mol.GetBondBetweenAtoms(atom_idx, neighbor_atom.GetIdx())
            for neighbor_atom in mapped_neighbors(mol.GetAtomWithIdx(atom_idx), as_indices=False)
    }
    
class AnnotatedReaction(ChemicalReaction):
    '''
    RDKit ChemicalReaction subclass with additional useful information about product atom and bond mappings and reaction naming
    Initialization must be done either via AnnotatedReaction.from_smarts, AnnotatedReaction.from_rdmols, or AnnotatedReaction.from_rxnfile
    '''
    # line number in .rxn file where (optional) name of reaction should be located (per the CTFile spec https://discover.3ds.com/sites/default/files/2020-08/biovia_ctfileformats_2020.pdf) 
    RXNNAME_LINE_NO : ClassVar[int] = 1
    RXNNAME_RE : ClassVar[re.Pattern] = re.compile(r'^\t*(?P<rxnname>.*?)\n$')

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)

    # INITIALIZATION/EXPORT TO SMARTS
    @classmethod
    def from_smarts(cls, rxn_smarts : str) -> 'AnnotatedReaction':
        '''Instantiate reaction from mapped SMARTS string'''
        return cls(ReactionFromSmarts(rxn_smarts.replace('#0', '*'))) # clean up any wild-atom conversion artifacts from porting a SMARTS through SMILES
    # NOTE : cannot analogous implement "from_smiles" classmethod, as rdChemreactions does not support initialization from SMILES (only SMARTS)

    def to_smarts(self) -> str:
        '''Export reaction as mapped SMARTS string''' # TODO : implement * -> R replacement here (rather than in rxn file I/O)
        return ReactionToSmarts(self).replace('#0', '*') # clean up any wild-atom conversion artifacts from porting a SMARTS through SMILES 

    @property
    def smarts(self) -> str:
        '''Mapped SMARTS string representation of reaction'''
        return self.to_smarts()

    def to_smiles(self) -> str:
        '''Export reaction as mapped SMILES string'''
        return ReactionToSmiles(self) # TODO : implement * -> R replacement here (rather than in rxn file I/O)

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

        return cls(ReactionFromMolecule(rxn_mol))

    # MDL RXN FILE I/O
    @classmethod
    @allow_string_paths
    def rxnname_from_rxnfile(cls, rxnfile_path : Union[str, Path]) -> str:
        '''Extract the reaction name from a (properly-formatted) .RXN file'''
        with rxnfile_path.open('r') as file:
            for i, line in enumerate(file):
                if i == cls.RXNNAME_LINE_NO:
                    return re.match(cls.RXNNAME_RE, line).group('rxnname')
                
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
        rxn = cls(ReactionFromRxnFile(rxnfile_path))
        rxn.rxnname = cls.rxnname_from_rxnfile(rxnfile_path)

        return rxn
    
    @allow_string_paths
    def to_rxnfile(self, rxnfile_path : Union[str, Path], wilds_to_R_groups : bool=True) -> None:
        '''Save reaction to an MDL .RXN file. Replaces ports with R-groups to enable proper loading'''
        rxn_block = ReactionToRxnBlock(self)
        if wilds_to_R_groups:
            rxn_block = rxn_block.replace('*', 'R')

        with rxnfile_path.open('w') as rxnfile:
            for i, line in enumerate(StringIO(rxn_block)):
                if (i == self.RXNNAME_LINE_NO):
                    rxnfile.write(f'\t\t\t{self.rxnname}\n')
                else:
                    rxnfile.write(line)

    # TRACING PROVENANCE OF CHEMICAL OBJECTS THROUGH THE REACTION
    ## ATOMS
    @cached_property
    def mapped_atom_info_by_map_number(self) -> dict[int, AtomTraceInfo]:
        '''Provenance info about the reactant origin and product destination of all mapped atoms, keyed by atom map numbers'''
        # 1) initialize and collect reactant side of atom info
        atom_info_aggregator : dict[int, dict[str, int]] = {}
        for reactant_idx, reactant_template in enumerate(self.GetReactants()):
            for reactant_atom in mapped_atoms(reactant_template, as_indices=False):
                map_number = reactant_atom.GetAtomMapNum()
                atom_info_aggregator[map_number] = {
                    'map_number' : map_number, # NOTE: need to make sure these match with the attribute names of AtomTraceInfo for now
                    'reactant_idx' : reactant_idx,
                    'reactant_atom_idx' : reactant_atom.GetIdx(),
                }
                
        # 2) fill product side info into existing data 
        for product_idx, product_template in enumerate(self.GetProducts()):
            for product_atom in mapped_atoms(product_template, as_indices=False):
                map_number = product_atom.GetAtomMapNum()
                atom_info_aggregator[map_number].update({ # DEVNOTE: don't need defaultdict here, since atom map numbers are bijective across reactions (by RDKit stipulation)
                    'product_idx' : product_idx,
                    'product_atom_idx' : product_atom.GetIdx(),
                })
              
        # 3) upnack collated info into (frozen) AtomTraceInfo objects
        return {
            map_number : AtomTraceInfo(**info_dict)
                for map_number, info_dict in atom_info_aggregator.items()
        }
    
    @cached_property
    def mapped_atom_info(self) -> set[AtomTraceInfo]:
        '''Compile provenance info about the reactant origin and product destination of all mapped atoms'''
        return set(self.mapped_atom_info_by_map_number.values())
    
    @cached_property
    def mapped_atom_info_by_product_idx(self) -> dict[str, set[AtomTraceInfo]]:
        '''
        Sets of provenance info about the reactant origin and product destination of all mapped atoms
        grouped by the indices of the products in which they appear
        '''
        atom_info_by_product_idx = defaultdict(set)
        for atom_info in self.mapped_atom_info:
            atom_info_by_product_idx[atom_info.product_idx].add(atom_info)

        return dict(atom_info_by_product_idx) # convert back to "regular" dict for typing

    @cached_property
    def reactive_atom_info(self) -> dict[int, AtomTraceInfo]:
        '''Compile reactant origin and product destination of all mapped atoms which are changed by the reaction'''
        reactive_atom_infos = {}
        for reactant_idx, reactant_atom_idxs in enumerate(self.GetReactingAtoms(mappedAtomsOnly=True)):
            reactant_template = self.GetReactantTemplate(reactant_idx)
            for reactant_atom_idx in reactant_atom_idxs:
                map_number = reactant_template.GetAtomWithIdx(reactant_atom_idx).GetAtomMapNum()
                reactive_atom_infos[map_number] = self.mapped_atom_info_by_map_number[map_number]
        
        return reactive_atom_infos
        
    ## BONDS
    @cached_property
    def mapped_bond_info(self) -> set[BondTraceInfo]:
        '''All provenance info on how bonds between mapped atoms (at least one of which is reactive) change over the reaction'''
        mapped_bond_infos = set()
        for map_number, atom_info in self.reactive_atom_info.items():
            reactant_neighbor_bonds = map_numbers_to_neighbor_bonds(self.GetReactantTemplate(atom_info.reactant_idx), atom_info.reactant_atom_idx)
            product_neighbor_bonds  = map_numbers_to_neighbor_bonds(self.GetProductTemplate(atom_info.product_idx), atom_info.product_atom_idx)
            
            # iterate over all mapped atoms that were, at any point in the reaction history, neighbors of the current reactive atom
            for neighbor_map_number in (reactant_neighbor_bonds.keys() | product_neighbor_bonds.keys()):
                # returns None either if the mapped atom is not present OR if it is present but not bonded to the target atom
                initial_bond : Optional[Bond] = reactant_neighbor_bonds.get(neighbor_map_number) 
                if initial_bond is None:
                    initial_bond_type = None
                else:
                    initial_bond_type = initial_bond.GetBondType()

                # returns None either if the mapped atom is not present OR if it is present but not bonded to the target atom
                final_bond : Optional[Bond] = product_neighbor_bonds.get(neighbor_map_number) 
                if final_bond is None:
                    final_bond_type = BondType.ZERO 
                    product_idx = None
                    product_bond_idx = None
                else:
                    final_bond_type = final_bond.GetBondType()
                    product_idx = atom_info.product_idx
                    product_bond_idx = final_bond.GetIdx()

                if initial_bond and final_bond:
                    bond_change_type = BondChange.UNCHANGED if (initial_bond_type == final_bond_type) else BondChange.MODIFIED
                elif initial_bond and not final_bond:
                    bond_change_type = BondChange.DELETED
                elif not initial_bond and final_bond:
                    bond_change_type = BondChange.ADDED

                mapped_bond_infos.add(
                    BondTraceInfo( 
                        map_nums=frozenset([map_number, neighbor_map_number]),
                        product_idx=product_idx,
                        product_bond_idx=product_bond_idx,
                        bond_change_type=bond_change_type,
                        initial_bond_type=initial_bond_type,
                        final_bond_type=final_bond_type,
                    )
                )
        return mapped_bond_infos

    @cached_property
    def mapped_bond_info_by_change_type(self) -> dict[Union[str, BondChange], set[BondTraceInfo]]:
        '''Mapped bond provenance information, grouped by the types of bond changes each bond experienced'''
        bond_info_by_change_type = defaultdict(set)
        for bond_info in self.mapped_bond_info:
            bond_info_by_change_type[bond_info.bond_change_type].add(bond_info)

        return dict(bond_info_by_change_type) # convert back to "regular" dict for typing

    @cached_property
    def mapped_bond_info_by_product_idx(self) -> dict[Union[str, BondChange], set[BondTraceInfo]]:
        '''Mapped bond provenance information, grouped by the index of the product the bond ends up in'''
        bond_info_by_product_idx = defaultdict(set)
        for bond_info in self.mapped_bond_info:
            bond_info_by_product_idx[bond_info.product_idx].add(bond_info)

        return dict(bond_info_by_product_idx) # convert back to "regular" dict for typing
    
    # REACTANT PATHFINDING
    def compile_functional_group_inventory(self, labeled_reactants : dict[L, Mol]) -> SymbolInventory[int, L]:
        '''
        Construct an inventory of numbers of functional groups (reactant templates) found in an ordered pool of reactant Mols,
        which can be evaluated as a DISCERNMENT-type problem to determine valid reactant ordering(s), if some exist
        
        Isomorphisms to the DISCERNMENT problem in this instance are as follows:
            "symbols"     <-> indices of functional groups, as-defined by reactant templates
            "target word" <-> the sequence of indices [n-1] = [0, 1, ..., n-1] where "n" is the number of reactant templates
            "magazine"    <-> an ordered collection of reactant molecules, which may contain any number of functional groups apiece each
            "word labels" <-> the output of the provided labeling method on the reactant Mol or, if no labeling method is provided, the index of the reactant in the sequence
        '''
        fn_group_sym_inv = defaultdict(Counter) # TODO: add more "native" mechanism to instantiate SymbolInventory directly from counts (rather than sequences of bins)
        for reactant_label, reactant in labeled_reactants.items(): # NOTE: placed first in case reactants are exhaustible (e.g. generator-like)
            for template_idx, reactant_template in enumerate(self.GetReactants()): # NOTE: can't simply replace this loop with a dict comprehension since any given key pair my not yet exist in the defaultdict
                fn_group_sym_inv[template_idx][reactant_label] = num_substruct_queries_distinct(reactant, reactant_template) # counts are the number of times a particular reactant template occurs in a monomer
        
        return SymbolInventory(fn_group_sym_inv)
        
    def enumerate_valid_reactant_orderings(
            self,
            reactant_pool : Union[Iterable[Mol], Mapping[L, Mol]],
            labeling_method : Optional[Callable[[Mol], L]]=None,
            as_mols : bool=True,
            allow_resampling : bool=False,
            deterministic : bool=True,
        ) -> Generator[Union[None, tuple[L], tuple[Mol]], None, None]:
        '''
        Enumerates all orderings of reactants compatible with the reactant templates defined in this reaction

        Yields:
        * a single NoneType if no such ordering exists
        * tuples of Chem.Mol objects if as_mols=True
        * tuples of molecule labels if as_mols=False. Labels can be:
            ** Indices of molecules in-sequence if labeling_method is not provided 
            ** Any molecule label is labeling_method is provided;
               Note that the labels provided must be valid keys for indexing reactant_pool to support as_mols
        
        If allow_resampling=False, each reactant will only be allowed to contribute exactly 1 of its functional groups 1 to any solution 
        '''
        # enforce dict-like reactant specification (with default implementations of label if none are given)
        if not isinstance(reactant_pool, Mapping):
            reactant_pool : dict[L, Mol] = {
                reactant_idx if (labeling_method is None) else labeling_method(reactant) : reactant # in absence of specialized labeling method, just use indices of molecules
                    for reactant_idx, reactant in enumerate(reactant_pool)
        }
        
        # shuffle order of reactants if some stochasticity is desired
        if not deterministic:
            reactant_labels = list(reactant_pool.keys())
            shuffle(reactant_labels)
            reactant_pool = {
                reactant_label : reactant_pool[reactant_label] 
                    for reactant_label in reactant_labels
            }

        # solve for all ordered reactant subsets compatible with the reaction 
        reactant_orderings_found : bool = False
        reactant_ordering_planner = DISCERNMENTSolver(
            self.compile_functional_group_inventory(reactant_pool)
        )
        for reactant_labels in reactant_ordering_planner.enumerate_choices(
            word=[templ_idx for templ_idx in range(self.GetNumReactantTemplates())], # use indices of reactants to allow direct lookup of reactants from solution
            unique_bins=not allow_resampling, # if one is OK with the same reactant being used for multiple functional groups, ignore its functional group multiplicities
        ):                                          
            yield tuple(reactant_labels) if not as_mols else tuple(reactant_pool[label] for label in reactant_labels) # by default, labels are reactant indices
            reactant_orderings_found = True # flag that at least one ordering is possible
        
        # if no solution exists, explicitly yield a single NoneType sentinel value (simplifies check for no solutions existing)
        if not reactant_orderings_found:  
            yield None 
    
    def valid_reactant_ordering(
            self,
            reactant_pool : Sequence[Mol],
            as_mols : bool=True,
            allow_resampling : bool=False,
            deterministic : bool=True,
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
            reactant_pool=reactant_pool,
            as_mols=as_mols,
            allow_resampling=allow_resampling,
            deterministic=deterministic,
            labeling_method=None, # DEVNOTE : for now, explicitly force reactant labeling to be index-based for this check to simplify downstream methods
        )) # DEVNOTE: enumeration guarantees a NoneType is returned when no solution exists (no edge-case handling needed!)
        
    def has_reactable_subset(self, reactant_pool : Sequence[Mol], allow_resampling : bool=False) -> bool:
        '''
        Determine if a sequence of reactants Mols contains any subset of Mols which are compatible with the reactant templates defined by this reaction
        If allow_resampling=False, each reactant will only be allowed to contribute exactly 1 of its functional groups 1 to any solution 
        '''
        return self.valid_reactant_ordering(
            reactant_pool=reactant_pool,
            as_mols=False,
            allow_resampling=allow_resampling,
            deterministic=True,
        ) is not None

    # RUNNING CHEMICAL REACTIONS
    def validate_reactants(self, reactants : Sequence[Mol], allow_resampling : bool=False) -> None:
        '''Check whether a collection of reactant Mols can be reacted with this reaction definition'''
        # can quickly discount a bad reactant sequence by a simple counting check, prior to the more expensive reactant order determination
        if (num_reactants_provided := len(reactants)) != (num_reactant_templates_required := self.GetNumReactantTemplates()):
            raise BadNumberReactants(f'{self.__class__.__name__} expected {num_reactant_templates_required} reactants, but {num_reactants_provided} were provided')
        
        # otherwise, look for any solutions to combinatorial reactant search
        if not self.has_reactable_subset(reactant_pool=reactants, allow_resampling=allow_resampling):  # check that the reactants are compatible with the reaction
            raise ReactantTemplateMismatch(f'Reactants provided to {self.__class__.__name__} are incompatible with reaction schema defined')
        
    def reactants_are_compatible(self, reactants : Sequence[Mol], allow_resampling : bool=False) -> bool:
        '''Determine whether a collection of reactants can be reacted with this reaction or not'''
        try:
            self.validate_reactants(reactants, allow_resampling=allow_resampling)
        except:
            return False
        else: 
            return True
    
    @staticmethod
    def apply_atom_info_to_product(
            product : Mol,
            product_atom_infos : Iterable['AtomTraceInfo'],
            reactants : Sequence[Mol],
            apply_map_labels : bool=True,
        ) -> None:
        '''Transfer props and (if requested) map number information from atoms in reactant Mols to their corresponding atoms in a product Mol
        Acts in-place on the "product" Mol instance'''
        for atom_info in product_atom_infos:
            product_atom = product.GetAtomWithIdx(atom_info.product_atom_idx)
            assert product_atom.HasProp('old_mapno') # precisely the mapped atoms in the reaction template will have this property set on the product 
            assert product_atom.GetIntProp('old_mapno') == atom_info.map_number # sanity check to make sure atom map identity has been preserved

            corresp_reactant_atom = reactants[atom_info.reactant_idx].GetAtomWithIdx(atom_info.reactant_atom_idx) # product_atom.GetIntProp('react_atom_idx'))
            copy_rdobj_props(from_rdobj=corresp_reactant_atom, to_rdobj=product_atom) # clone props from corresponding atom in reactant
            
            if apply_map_labels:
                product_atom.SetAtomMapNum(atom_info.map_number)
                product_atom.ClearProp('old_mapno')
                
    @staticmethod
    def apply_bond_info_to_product(
            product : Mol,
            product_bond_infos : Iterable['BondTraceInfo'],
        ) -> None:
        '''Mark any changed bonds with bond props and clean up bond type info in places where bonds get modified
        Acts in-place on the "product" Mol instance'''
        for bond_info in product_bond_infos:
            product_bond = product.GetBondWithIdx(bond_info.product_bond_idx)
            
            if bond_info.bond_change_type == BondChange.ADDED: # explicictly labl new or changed bonds
                product_bond.SetProp(BOND_CHANGE_PROPNAME, BondChange.ADDED)
    
            if bond_info.bond_change_type == BondChange.MODIFIED:
                product_bond.SetProp(BOND_CHANGE_PROPNAME, BondChange.MODIFIED)
                assert(product_bond.GetBeginAtom().HasProp('_ReactionDegreeChanged')) # double check reaction agrees the bond has changed
                assert(product_bond.GetEndAtom().HasProp('_ReactionDegreeChanged')) 

                if bond_info.final_bond_type not in (BondType.ZERO, BondType.UNSPECIFIED): # if an explicit bond type is defined in the template...
                    product_bond.SetBondType(bond_info.final_bond_type) # ...set the product's bond type to what it *should* be from the reaction schema

    @sanitizable_mol_outputs
    def react(
            self,
            reactants : Sequence[Mol],
            repetitions : int=1,
            keep_map_labels : bool=True,
            _suppress_reactant_validation : bool=False,
        ) -> Generator[Mol, None, None]:
        '''
        Execute reaction over a collection of reactants and generate product molecule(s)
        Does not require reactants to match the ORDER of the expected reactant templates by default 
        (only to have the correct number of reactants)
        '''
        # validate reactants
        if not _suppress_reactant_validation:
            self.validate_reactants(reactants, allow_resampling=False) # NOTE: reactants here are expected to match templates 1:1, so resampling should never be necessary
        
        # validate and label reactants
        reactants = [Mol(reactant) for reactant in reactants] # make a copy to avoid preserve read-onlyness of inputted reactant Mols
        for reactant_idx, reactant in enumerate(reactants):
            for atom in reactant.GetAtoms():
                atom.SetIntProp(REACTANT_INDEX_PROPNAME, reactant_idx) # label reactant atoms with their respective reactant IDs 
        
        # iterate over raw RDKit products, sanitizing and injecting information before yielding
        raw_products = self.RunReactants(reactants, maxProducts=repetitions) # obtain unfiltered RDKit reaction output. TODO : generalize to work when more than 1 repetition is requested
        for product_idx, product in enumerate(chain.from_iterable(raw_products)): # clean up products into a usable form
            self.apply_atom_info_to_product( # copy reactant atom props over to product atoms
                product,
                product_atom_infos=self.mapped_atom_info_by_product_idx[product_idx],
                reactants=reactants,
                apply_map_labels=keep_map_labels,
            )
            self.apply_bond_info_to_product( # indicate additions or modifications to product bonds
                product,
                product_bond_infos=self.mapped_bond_info_by_product_idx[product_idx],
            )
            yield product