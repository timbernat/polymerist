'''For fragmenting molecules by reaction and residue information'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Generator
from abc import ABC, abstractmethod
from itertools import combinations

from rdkit import Chem
from rdkit.Chem import rdqueries, Mol

from .reactions import RxnProductInfo
from ..labeling import bondwise


# ABSTRACT BASE FOR FRAGMENTATION STRATEGIES
class IntermonomerBondIdentificationStrategy(ABC):
    '''Abstract base for Intermonomer Bond Identification Strategies for fragmentation during in-silico polymerization'''
    @abstractmethod
    def locate_intermonomer_bonds(self, product : Mol, product_info : RxnProductInfo) -> Generator[int, None, None]:
        '''Generates the indices of all identified inter-monomer bonds by molecule'''
        pass

    def produce_fragments(self, product : Mol, product_info : RxnProductInfo, separate : bool=True):
        '''Apply break all bonds identified by this IBIS algorithm and return the resulting fragments'''
        fragments = Chem.FragmentOnBonds(
            mol=product,
            bondIndices=self.locate_intermonomer_bonds(product, product_info) # TODO : check that the multiplicity of any bond to cut is no greater than the bond order
        ) # TODO : add config for "dummyLabels" arg to support port flavor setting
        if separate:
            return Chem.GetMolFrags(fragments, asMols=True, sanitizeFrags=False) # avoid disruptive sanitization (can be done in post)
        return fragments # if separation is not requested, return as single unfragmented molecule object
IBIS = IntermonomerBondIdentificationStrategy # shorthand alias for convenience

# CONCRETE IMPLEMENTATIONS
DUMMY_PROP_QUERY = rdqueries.HasPropQueryAtom('was_dummy') # heavy atom which was converted from a dummy atom in a reaction
HEAVY_FORMER_LINKER_QUERY = Chem.MolFromSmarts('A')
HEAVY_FORMER_LINKER_QUERY.GetAtomWithIdx(0).ExpandQuery(DUMMY_PROP_QUERY) # cast as Mol to allow for quick check via GetSubstructMatch

class ReseparateRGroups(IBIS):
    '''IBIS which cleaves any new bonds formed between atoms that were formerly the start of an R-group in the reaction template'''
    def locate_intermonomer_bonds(self, product: Mol, product_info : RxnProductInfo) -> Generator[int, None, None]:
        possible_bridgehead_ids = [atom_id for match in product.GetSubstructMatches(HEAVY_FORMER_LINKER_QUERY) for atom_id in match]
        for new_bond_id in product_info.new_bond_ids_to_map_nums.keys():                     # for each newly formed bond...
            for bridgehead_id_pair in combinations(possible_bridgehead_ids, 2):                   # ...find the most direct path between bridgehead atoms...
                if new_bond_id in bondwise.get_shortest_path_bonds(product, *bridgehead_id_pair): # ...and check if the new bond lies along it
                    yield new_bond_id

class ReseparateRGroupsUnique(IBIS):
    '''IBIS which cleaves any new bonds formed between atoms that were formerly the start of an R-group in the reaction template exactly once each'''
    def locate_intermonomer_bonds(self, product: Mol, product_info : RxnProductInfo) -> Generator[int, None, None]:
        possible_bridgehead_ids : list[int] = [atom_id for match in product.GetSubstructMatches(HEAVY_FORMER_LINKER_QUERY) for atom_id in match]
        bonds_already_cut : set[int] = set()

        for bridgehead_id_pair in combinations(possible_bridgehead_ids, 2):                                     # for every pair of R-group bridgehead atoms...
            shortest_path_bond_ids : list[int] = bondwise.get_shortest_path_bonds(product, *bridgehead_id_pair) # find the shortest bond path between the bridgeheads... 
            for new_bond_id in product_info.new_bond_ids_to_map_nums.keys():                                    # then for each newly formed bond...
                if (new_bond_id in shortest_path_bond_ids) and (new_bond_id not in bonds_already_cut):          # if the new bond lies along the path and has not already been selected for cutting...
                    yield new_bond_id                                                                           # select the new bond for cutting...
                    bonds_already_cut.add(new_bond_id)                                                          # and mark it as visited to avoid duplicate cuts