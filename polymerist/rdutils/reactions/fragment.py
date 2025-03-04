'''For fragmenting molecules by reaction and residue information'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Callable, ClassVar, Generator, Optional, ParamSpec
P = ParamSpec('P')

from abc import ABC, abstractmethod
from itertools import combinations

from rdkit import Chem
from rdkit.Chem import rdqueries, Mol, Bond
from networkx import bridges

from .reactinfo import BOND_CHANGE_PROPNAME, BondChange
from ..selection import (
    atom_adjoins_linker,
    BondCondition,
    bonds_by_condition,
    bond_condition_by_atom_condition_factory,
    logical_or,
)
from ..rdgraphs import chemical_graph
from ...genutils.iteration import sliding_window


# HELPER FUNCTIONS
## BRIDGEHEAD ATOM PATHFINDING
HEAVY_FORMER_LINKER_QUERY_ATOM : Chem.QueryAtom = rdqueries.HasPropQueryAtom('was_dummy') # query for atoms with the property "was_dummy" set
HEAVY_FORMER_LINKER_QUERY_ATOM.ExpandQuery(rdqueries.AAtomQueryAtom()) # expand query to only match heavy atoms
def bridgehead_atom_ids(product : Chem.Mol) -> Generator[int, None, None]:
    '''
    Generates the indices of all atoms in a reaction product which were tagged
    as R-group bridgehead (i.e. wild) atoms in the reaction template definition
    '''
    for bh_atom in product.GetAtomsMatchingQuery(HEAVY_FORMER_LINKER_QUERY_ATOM):
        yield bh_atom.GetIdx()
        
def get_shortest_path_bonds(rdmol : Mol, start_atom_idx : int, end_atom_idx : int) -> list[int]:
    '''Returns bond indices along shortest path between two atoms in a Mol'''
    return [
        rdmol.GetBondBetweenAtoms(*atom_id_pair).GetIdx()
            for atom_id_pair in sliding_window(Chem.GetShortestPath(rdmol, start_atom_idx, end_atom_idx), n=2)
    ]
    
## BOND CUTTING CRITERIA
def bond_is_newly_formed(bond : Bond) -> bool:
    '''Bond condition checking if a bond was newly formed in a reaction
    (i.e. present in the product(s) but not the reactant(s))'''
    return bond.HasProp(BOND_CHANGE_PROPNAME) and (bond.GetProp(BOND_CHANGE_PROPNAME) == BondChange.ADDED) 

bond_adjoins_linker : BondCondition = bond_condition_by_atom_condition_factory(
    atom_condition=atom_adjoins_linker,
    binary_operator=logical_or,
)

# FRAGMENTATION STRATEGIES
class IntermonomerBondIdentificationStrategy(ABC):
    '''Abstract base for Intermonomer Bond Identification Strategies for fragmentation during in-silico polymerization'''
    @abstractmethod
    def _locate_intermonomer_bonds(self, product : Mol) -> Generator[int, None, None]:
        '''
        Generates the indices of all identified inter-monomer bonds by molecule
        MUST BE IMPLEMENTED in order to define behavior of fragmentation strategy
        '''
        pass

    def locate_intermonomer_bonds(self, product : Mol) -> Generator[int, None, None]:
        '''Generates the indices of all identified inter-monomer bonds by molecule, no more than once each'''
        bonds_already_cut : set[int] = set()
        for bond_id in self._locate_intermonomer_bonds(product):
            if bond_id not in bonds_already_cut: # bond cleavage must be idempotent, to avoid attempting to cut bonds which no longer exist
                yield bond_id
                bonds_already_cut.add(bond_id)   # mark bond as visited to avoid duplicate cuts

    def produce_fragments(self, product : Mol, separate : bool=True):
        '''Apply break all bonds identified by this IBIS algorithm and return the resulting fragments'''
        fragments = Chem.FragmentOnBonds(
            mol=product,
            bondIndices=self.locate_intermonomer_bonds(product) # TODO : check that the multiplicity of any bond to cut is no greater than the bond order
        ) # TODO : add config for "dummyLabels" arg to support port flavor setting
        if separate:
            return Chem.GetMolFrags(fragments, asMols=True, sanitizeFrags=False) # avoid sanitizing molecules prematurely (can be done in post)
        return fragments # if separation is not requested, return as single unfragmented molecule object
IBIS = IntermonomerBondIdentificationStrategy # shorthand alias for convenience
   
## CONCRETE IMPLEMENTATIONS
class ReseparateRGroups(IntermonomerBondIdentificationStrategy):
    '''IBIS which cleaves any new bonds formed between atoms that were formerly the start of an R-group in the reaction template'''
    def _locate_intermonomer_bonds(self, product: Mol) -> Generator[int, None, None]:
        new_bond_ids : set[int] = set(bonds_by_condition(product, condition=bond_is_newly_formed, as_indices=True, as_pairs=False, negate=False))
        for bridgehead_id_pair in combinations(bridgehead_atom_ids(product), 2):         # for every pair of R-group bridgehead atoms...
            for new_bond_id in new_bond_ids:                                             # find the path(s) with fewest bonds between the bridgeheads... 
                if new_bond_id in get_shortest_path_bonds(product, *bridgehead_id_pair): # and select for cutting any newly-formed bonds found along that path
                    yield new_bond_id
                    
class CutMinimumWeightBondsStrategy(IntermonomerBondIdentificationStrategy):
    '''Subtype of IBIS which chooses bonds as solutions to minimum graph cut problem
    All bonds in a molecule are given some base cost to cut, then a discount is applied to some bonds by specified conditions
    
    Cuts are then made (no more than once) on the lowest cost bond(s) separating each pair of R-groups'''
    _DEFAULT_BASE_BOND_COST : ClassVar[float] = 1.0
    _DEFAULT_BOND_DISCOUNTS : ClassVar[dict[str, tuple[int, Callable[[Mol], tuple[int, int]]]]] = {
        'BRIDGE_BOND'       : (0.5  , lambda mol : bridges(chemical_graph(mol))), # networkx kindly already returns these as node index pairs
        'NEW BOND'          : (0.25 , lambda mol : bonds_by_condition(mol, bond_is_newly_formed, as_indices=True, as_pairs=True, negate=False)),
        'RGROUP-FREE BOND'  : (0.125, lambda mol : bonds_by_condition(mol, bond_adjoins_linker , as_indices=True, as_pairs=True, negate=True)),
    }
    _DEFAULT_COST_KEYWORD : ClassVar[str] = 'capacity'
    
    def __init__(
            self,
            base_bond_cost : Optional[float]=None,
            bond_discounts : Optional[dict[str, tuple[int, Callable[[Mol], tuple[int, int]]]]]=None,
            bond_cost_keyword : Optional[str]=None,
            *args : P.args,
            **kwargs : P.kwargs,
        ):
        super().__init__(*args, **kwargs)
        
        # set default values as necessary
        if base_bond_cost is None:
            base_bond_cost = self._DEFAULT_BASE_BOND_COST
            
        if bond_discounts is None:
            bond_discounts = self._DEFAULT_BOND_DISCOUNTS
            
        # validate bond weights provided
        if base_bond_cost <= 0.0:
            raise ValueError(f'Must provide positive base bond cutting cost (provided {base_bond_cost})')
        
        max_discount = sum(discount for discount, _ in bond_discounts.values())
        if max_discount > base_bond_cost:
            raise ValueError(f'Maximum possible discount of {max_discount} cannot exceed base cost {base_bond_cost} for cutting a bond')
        
        # set attributes
        self.base_bond_cost = base_bond_cost
        self.bond_discounts = bond_discounts
        self.bond_cost_keyword = bond_cost_keyword
        
    def _locate_intermonomer_bonds(self, product: Mol) -> Generator[int, None, None]:
        raise NotImplemented
                                                    