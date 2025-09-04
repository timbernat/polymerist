'''Tools for simplifying the construction of reaction templates'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from dataclasses import dataclass, field
from typing import Iterable, Optional, Union

from rdkit.Chem.rdchem import Mol, RWMol
from rdkit.Chem.rdmolops import GetMolFrags, SanitizeMol, SanitizeFlags, SANITIZE_ALL

from .reactions import AnnotatedReaction
from ..chemlabel import relabel_map_nums

from ..bonding._bonding import combined_rdmol
from ..bonding.permutation import swap_bonds


@dataclass # TODO : make JSONifiable
class ReactionAssembler:
    '''Class for producing reaction templates from reactants, bond derangements'''
    reactive_groups  : Iterable[Mol]
    byproducts       : Iterable[Mol]              = field(default_factory=list)
    bond_derangement : dict[int, tuple[int, int]] = field(default_factory=dict)
    rxn_name         : str = ''
    # TODO : incorporate stereo

    @property
    def reactants(self) -> RWMol:
        '''Combine cached reactive groups into single, contiguously-numbered Mol for manipulation'''
        return combined_rdmol(
            *self.reactive_groups,
            # NOTE: unique mapping and editablility are essential for well-defined behavior throughout!
            assign_map_nums=True,
            editable=True,
        )

    def products(self,
            show_steps : bool=False,
            sanitize_ops : SanitizeFlags=SANITIZE_ALL,
            bond_breakage_marker : str='--x->',
            bond_formation_marker : str='---->',
        ) -> Mol:
        '''Generate the product template defined by the provided reactants and bond derangement'''
        if not self.bond_derangement:
            raise ValueError('Must provide non-empty bond derangement')

        # 2) defining and swapping bonds to form product
        products = swap_bonds(
            self.reactants, # NOTE: this should be an RWMol (not just Mol) instance!
            bond_derangement=self.bond_derangement,
            show_steps=show_steps,
            bond_breakage_marker=bond_breakage_marker,
            bond_formation_marker=bond_formation_marker
        ) # create editable Mol
        SanitizeMol(products, sanitizeOps=sanitize_ops)

        return products
    
    def products_by_importance(
            self,
            combined : bool=True,
            show_steps : bool=False,
            sanitize_ops : SanitizeFlags=SANITIZE_ALL,
        ) -> tuple[Union[Optional[Mol], list[Mol]], Union[Optional[Mol], list[Mol]]]:
        '''Partition reaction products into major and minor/byproducts, each returned as a single Combined Mol'''
        product_partition = main_products, byproducts = [], []
        for product in GetMolFrags(self.products(show_steps=show_steps, sanitize_ops=sanitize_ops), asMols=True):
            for side_query_mol in self.byproducts:
                if product.HasSubstructMatch(side_query_mol) and (product.GetNumAtoms() == side_query_mol.GetNumAtoms()):
                    byproducts.append(product)
                    break # assumes uniquely-provided byproducts (a reasonable assumption)
            else:
                main_products.append(product)

        if not combined:
            return product_partition
        
        return [ # implicit else
            combined_rdmol(*mol_list, assign_map_nums=False, editable=False) if mol_list else None
                for mol_list in product_partition
        ]
    
    @property # TODO : add assertion that reactant and product map num sets match
    def byproduct_map_nums(self) -> tuple[set, set]:
        '''Partitions map numbers present in the product(s) by whether or not they belong to a collection of side products
        Returns a set of map numbers NOT in a side product and set set which are'''
        return [
            set(atom.GetAtomMapNum() for atom in product_mol.GetAtoms()) if (product_mol is not None) else set()
                for product_mol in self.products_by_importance(combined=True, show_steps=False, sanitize_ops=SANITIZE_ALL) # CRITICAL that mols be combined here
        ]

    @property
    def byproduct_relabeling(self) -> dict[int, int]:
        '''Determine a relabeling of the N map numbers in the reaction template which are not any of the given side products
        Relabeling is in standard form, i.e. new labels are taken from the first N natural numbers with the order of labels preserved'''
        map_nums_to_keep, map_nums_to_clear = self.byproduct_map_nums
        
        relabeling = {}
        for i, map_num in enumerate(map_nums_to_keep):
            relabeling[map_num] = (i + 1) # generate standard labeling of preserved map nums

        for map_num in map_nums_to_clear:
            relabeling[map_num] = 0 # unset labels on the unkept atoms

        return relabeling

    def assemble_rxn(self, show_steps : bool=False, sanitize_ops : SanitizeFlags=SANITIZE_ALL,) -> AnnotatedReaction:
        '''Assemble MDL rxn template from information stored in self'''
        reactants = self.reactants
        products, byproducts = self.products_by_importance(combined=True, show_steps=show_steps, sanitize_ops=sanitize_ops)

        if byproducts is not None:
            relabeling = self.byproduct_relabeling
            reactants = relabel_map_nums(reactants, relabeling=relabeling, in_place=False)
            products  = relabel_map_nums(products , relabeling=relabeling, in_place=False)

        rxn = AnnotatedReaction.from_rdmols(reactant_templates=[reactants], product_templates=[products])
        rxn.Initialize()
        rxn.rxnname = self.rxn_name
        
        num_warnings, num_errors = rxn.Validate()
        if num_errors != 0:
            raise ValueError('Issues with reaction definition')

        return rxn