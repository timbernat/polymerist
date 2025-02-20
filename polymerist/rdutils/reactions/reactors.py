'''Classes for implementing reactions with respect to some set of reactant RDKit Mols'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Generator, Iterable
from dataclasses import dataclass
from itertools import chain

from rdkit import Chem
from rdkit.Chem.rdchem import Mol, BondType
from rdkit.Chem.rdmolops import SanitizeFlags, SANITIZE_ALL

from .reactexc import BadNumberReactants, ReactantTemplateMismatch
from .reactions import AnnotatedReaction, BondChange, REACTANT_INDEX_PROPNAME, BOND_CHANGE_PROPNAME
from .fragment import IBIS, ReseparateRGroups

from ..rdprops import copy_rdobj_props
from ..rdprops.atomprops import clear_atom_props
from ..rdprops.bondprops import clear_bond_props
from ..labeling.molwise import clear_atom_map_nums


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

    # POST-REACTION CLEANUP METHODS

    # REACTION EXECUTION METHODS
    def react(
            self,
            reactants : Iterable[Mol],
            repetitions : int=1,
            keep_map_labels : bool=True,
            sanitize_ops : SanitizeFlags=SANITIZE_ALL,
        ) -> list[Mol]:
        '''Execute reaction over a collection of reactants and generate product molecule(s)
        Does NOT require the reactants to match the order of the reaction template (only that some order fits)'''
        # can quickly discount a bad reactant sequence by a simple counting check, prior to the more expensive reactant order determination
        if (num_reactants_provided := len(reactants)) != (num_reactant_templates_required := self.rxn_schema.GetNumReactantTemplates()):
            raise BadNumberReactants(f'{self.__class__.__name__} expected {num_reactant_templates_required} reactants, but {num_reactants_provided} were provided')
        
        reactants = self.rxn_schema.valid_reactant_ordering(reactants, as_mols=True) # check that the reactants are compatible with the reaction
        if reactants is None:
            raise ReactantTemplateMismatch(f'Reactants provided to {self.__class__.__name__} are incompatible with reaction schema defined')
        
        # label reactant atoms with their respective reactant IDs 
        reactants = [Chem.Mol(reactant) for reactant in reactants] # make a copy to avoid preserve read-onliness of inputted reactant Mols
        for reactant_idx, reactant in enumerate(reactants):
            for atom in reactant.GetAtoms():
                atom.SetIntProp(REACTANT_INDEX_PROPNAME, reactant_idx)
        
        # iterate over raw RDKit products, sanitizing and injecting information before yielding
        products : list[Mol] = []
        raw_products = self.rxn_schema.RunReactants(reactants, maxProducts=repetitions) # obtain unfiltered RDKit reaction output. TODO : generalize to work when more than 1 repetition is requested
        for product_idx, product in enumerate(chain.from_iterable(raw_products)): # clean up products into a usable form
            ## copy reactant atom props over to product atoms
            for product_atom in product.GetAtoms():
                if not product_atom.HasProp('old_mapno'): # precisely the mapped atoms in the reaction template will have this property set on the product 
                    continue
                
                map_number = product_atom.GetIntProp('old_mapno')
                atom_info = self.rxn_schema.mapped_atom_info_by_map_number[map_number]
                copy_rdobj_props( # clone props from corresponding atom in reactant
                    from_rdobj=reactants[atom_info.reactant_idx].GetAtomWithIdx(atom_info.reactant_atom_idx),
                    # from_rdobj=reactants[atom_info.reactant_idx].GetAtomWithIdx(product_atom.GetIntProp('react_atom_idx')),
                    to_rdobj=product_atom
                )
                if keep_map_labels:
                    product_atom.SetAtomMapNum(map_number)
                    product_atom.ClearProp('old_mapno')

            ## label added and modified bonds + clean up any modified bonds
            for bond_info in self.rxn_schema.mapped_bond_info_by_product_idx[product_idx]:
                assert bond_info.product_idx == product_idx # good sanity check, also helps to catch NoneType values for DELETED bonds
                product_bond = product.GetBondWithIdx(bond_info.product_bond_idx)
                bond_change_type : BondChange = bond_info.bond_change_type
                
                if bond_change_type in (BondChange.ADDED, BondChange.MODIFIED): # explicictly labl new or changed bonds
                    product_bond.SetProp(BOND_CHANGE_PROPNAME, bond_change_type)
     
                if bond_change_type == BondChange.MODIFIED:
                    ### double check that the reaction agrees that the bond has changed
                    assert(product_bond.GetBeginAtom().HasProp('_ReactionDegreeChanged')) 
                    assert(product_bond.GetEndAtom().HasProp('_ReactionDegreeChanged')) 

                    if bond_info.final_bond_type not in (BondType.ZERO, BondType.UNSPECIFIED): # if an explicit bond type is defined in the template...
                        product_bond.SetBondType(bond_info.final_bond_type) # ...set the product's bond type to what it *should* be from the reaction schema
            
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
                adducts = self.react(reactants, repetitions=1, sanitize_ops=sanitize_ops) # can't clear properties yet, otherwise intermonomer bond finder would have nothing to work with
            except ReactantTemplateMismatch:
                break
            
            fragments : list[Mol] = []
            for product in adducts: # DEVNOTE: consider doing fragmentation on the combined molecule made up of all products?
                for fragment in fragment_strategy.produce_fragments(product, separate=True):
                    clear_atom_props(fragment, in_place=True) # essential to avoid reaction mapping info from prior steps from contaminating future ones
                    clear_bond_props(fragment, in_place=True)
                    fragments.append(fragment)

                if clear_map_nums: # NOTE : CRITICAL that this be done after fragmentation step, which RELIES on map numbers being present
                    clear_atom_map_nums(product, in_place=True)
                    
            yield adducts, fragments # yield the adduct Mol and any subsequent resulting reactive fragments
            reactants = fragments # set fragments from current round of polymerization as reactants for next round
            
    def propagate_iterative(self) -> None:
        ...