'''Classes for implementing reactions with respect to some set of reactant RDKit Mols'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Generator, Iterable
from dataclasses import dataclass
from itertools import chain

from rdkit import Chem
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdmolops import SanitizeFlags, SANITIZE_ALL

from .reactexc import BadNumberReactants, ReactantTemplateMismatch
from .reactions import AnnotatedReaction
from .fragment import IBIS, ReseparateRGroups

from ..rdprops.atomprops import atom_ids_with_prop, clear_atom_props
from ..rdprops.bondprops import clear_bond_props
from ..labeling.molwise import clear_atom_map_nums, get_bond_by_map_num_pair


# REACTOR BASE CLASS 
@dataclass
class Reactor:
    '''Class for executing a reaction template on collections of RDKit Mol "reactants"'''
    rxn_schema : AnnotatedReaction

    # PRE-REACTION PREPARATION METHODS
    def _activate_reaction(self) -> None:
        '''Check that the reaction schema provided is well defined and initialized'''
        pass

    def __post_init__(self) -> None:
        '''Pre-processing of reaction and reactant Mols'''
        self._activate_reaction()

    @staticmethod
    def _label_reactants(rxn : AnnotatedReaction, reactants : Iterable[Mol]) -> None:
        '''Assigns "reactant_idx" Prop to all reactants to help track where atoms go during the reaction'''
        for i, reactant in enumerate(reactants):
            for atom in reactant.GetAtoms():
                atom.SetIntProp(rxn._atom_ridx_prop_name, i)

    # POST-REACTION CLEANUP METHODS
    @staticmethod
    def _relabel_reacted_atoms(rxn : AnnotatedReaction, product : Mol) -> None:
        '''Re-assigns "reactant_idx" Prop to modified reacted atoms to re-complete atom-to-reactant numbering'''
        for atom in product.GetAtoms():
            if atom.HasProp('old_mapno'):
                map_num = atom.GetIntProp('old_mapno')
                atom.SetIntProp(rxn._atom_ridx_prop_name, rxn.mapped_atom_info[map_num].reactant_idx)
                atom.SetAtomMapNum(map_num)

    @staticmethod
    def _clean_up_bond_orders(rxn : AnnotatedReaction, product : Mol, product_idx : int) -> None:
        '''Ensure bond order changes specified by the reaction are honored by RDKit'''
        product_info = rxn.product_info_maps[product_idx]
        product_template = rxn.GetProductTemplate(product_idx)
        
        for prod_templ_bond_id, map_num_pair in product_info.mod_bond_ids_to_map_nums.items():
            target_bond = product_template.GetBondWithIdx(prod_templ_bond_id)
            product_bond = get_bond_by_map_num_pair(product, map_num_pair, as_bond=True)

            assert(product_bond.GetBeginAtom().HasProp('_ReactionDegreeChanged')) 
            assert(product_bond.GetEndAtom().HasProp('_ReactionDegreeChanged')) # double check that the reaction agrees that the bond has changed

            product_bond.SetBondType(target_bond.GetBondType()) # set bond type to what it *should* be from the reaction schema

    @staticmethod
    def _label_new_and_modified_bonds(rxn : AnnotatedReaction, product : Mol, product_idx : int) -> None:
        '''Mark any bonds which were changed or added in the product'''
        product_info = rxn.product_info_maps[product_idx]
        
        _bond_change_labellers : dict[str, dict[int, tuple[int, int]]] = {
            'new_bond' : product_info.new_bond_ids_to_map_nums,
            'modified_bond' : product_info.mod_bond_ids_to_map_nums,
        }
        for bond_prop_value, changed_bond_dict in _bond_change_labellers.items():
            for map_num_pair in changed_bond_dict.values():
                product_bond = get_bond_by_map_num_pair(product, map_num_pair, as_bond=True)
                product_bond.SetProp(rxn._bond_change_prop_name, bond_prop_value)

    # REACTION EXECUTION METHODS
    def react(self, reactants : Iterable[Mol], repetitions : int=1, clear_props : bool=False, sanitize_ops : SanitizeFlags=SANITIZE_ALL) -> list[Mol]:
        '''Execute reaction over a collection of reactants and generate product molecule(s)
        Does NOT require the reactants to match the order of the reaction template (only that some order fits)'''
        # can quickly discount a bad reactant sequence by a simple counting check, prior to the more expensive reactant order determination
        if (num_reactants_provided := len(reactants)) != (num_reactant_templates_required := self.rxn_schema.GetNumReactantTemplates()):
            raise BadNumberReactants(f'{self.__class__.__name__} expected {num_reactant_templates_required} reactants, but {num_reactants_provided} were provided')
        
        reactants = self.rxn_schema.valid_reactant_ordering(reactants, as_mols=True) # check that the reactants are compatible with the reaction
        if reactants is None:
            raise ReactantTemplateMismatch(f'Reactants provided to {self.__class__.__name__} are incompatible with reaction schema defined')
        
        reactants = [Chem.Mol(reactant) for reactant in reactants] # make a copy to avoid preserve read-onliness of inputted reactant Mols
        self._label_reactants(
            self.rxn_schema,
            reactants=reactants,
        ) # assign reactant indices (not in-place)
        
        products : list[Mol] = []
        raw_products = self.rxn_schema.RunReactants(reactants, maxProducts=repetitions) # obtain unfiltered RDKit reaction output. TODO : generalize to work when more than 1 repetition is requested
        for i, product in enumerate(chain.from_iterable(raw_products)): # clean up products into a usable form
            self._relabel_reacted_atoms(
                self.rxn_schema,
                product=product,
            )
            self._clean_up_bond_orders(
                self.rxn_schema,
                product=product,
                product_idx=i,
            )
            self._label_new_and_modified_bonds(
                self.rxn_schema,
                product=product,
                product_idx=i,
            )
            if clear_props:
                clear_atom_props(product, in_place=True)
            Chem.SanitizeMol(product, sanitizeOps=sanitize_ops) # perform sanitization as-specified by the user
                
            products.append(product)
        return products


# REACTOR SUBCLASSES
@dataclass
class PolymerizationReactor(Reactor):
    '''Reactor which exhaustively generates monomers fragments according to a given a polymerization mechanism'''
    def propagate(
        self,
        monomers : Iterable[Mol],
        fragment_strategy : IBIS=ReseparateRGroups(),
        clear_map_nums : bool=True,
        sanitize_ops : SanitizeFlags=SANITIZE_ALL,
     ) -> Generator[tuple[list[Mol], list[Mol]], None, None]:
        '''Keep reacting and fragmenting a pair of monomers until all reactive sites have been reacted
        Returns fragment pairs at each step of the chain propagation process'''
        reactants = monomers # initialize reactive pair with monomers
        while True: # check if the reactants can be applied under the reaction template
            try:
                adducts = self.react(reactants, repetitions=1, clear_props=False, sanitize_ops=sanitize_ops) # can't clear properties yet, otherwise intermonomer bond finder would have nothing to work with
            except ReactantTemplateMismatch:
                break
            
            fragments : list[Mol] = []
            for product in adducts: # DEVNOTE: consider doing fragmentation on the combined molecule made up of all products?
                fragments.extend( # list extension preserves insertion order at each step
                    clear_bond_props(clear_atom_props(fragment, in_place=False), in_place=False) # essential to avoid reaction mapping info from prior steps from contaminating future ones
                        for fragment in fragment_strategy.produce_fragments(
                            product,
                            separate=True
                        )
                )
                # NOTE : CRITICAL that this be done after fragmentation step, which RELIES on map numbes being present
                if clear_map_nums:
                    clear_atom_map_nums(product, in_place=True)
            yield adducts, fragments # yield the adduct Mol and any subsequent resulting reactive fragments
            reactants = fragments # set fragments from current round of polymerization as reactants for next round