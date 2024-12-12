'''Tools for generating and manipulating monomer connectivity graphs'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Any, ClassVar, Generator, Optional, Sequence, Union

import networkx as nx
from itertools import product as cartesian_product

from ..genutils.iteration import asiterable
from ..genutils.sequences.discernment import DISCERNMENTSolver
from ..genutils.textual.delimiters import validate_braces
from ..genutils.fileutils.jsonio.serialize import TypeSerializer


class MonomerGraph(nx.Graph):
    '''A graph representation of the connectivity of monomer fragments in a polymer topology'''
    MONOMER_NAME_ATTR : ClassVar[str]            = 'monomer_name'     # node attribute name assigned to monomer names
    FLAVOR_DICT_ATTR  : ClassVar[dict[int, int]] = 'neighbor_flavors' # node attribute name assigned to outgoing flavors for bonds to neighbor ports
    BONDTYPE_ATTR     : ClassVar[str]            = 'bondtype'         # edge attribute name assigned to bond type annotations


    # node and edge attributes
    def get_monomer_name_at_node_index(self, node_idx : int) -> Optional[str]:
        '''Recover the assigned monomer name for the node at the given index'''
        return self.nodes[node_idx].get(self.MONOMER_NAME_ATTR)
    monomer_name = get_monomer_name = get_monomer_name_at_node_idx = get_monomer_name_at_node_index

    def get_flavor_dict_at_node_index(self, node_idx : int) -> Optional[dict[int, int]]:
        '''Recover the assigned dictionary of neighbor flavors for the node at the given index'''
        return self.nodes[node_idx].get(self.FLAVOR_DICT_ATTR)
    flavor_dict = get_flavor_dict = get_flavor_dict_at_node_idx = get_flavor_dict_at_node_index


    # connectivity properties
    @property
    def num_monomers(self) -> int:
        '''Number of monomer units represented in the current polymer'''
        return self.number_of_nodes()

    @property
    def is_unbranched(self) -> bool:
        '''Whether the monomer graph represents straight chain(s) without branching'''
        return all(node_deg <= 2 for node_id, node_deg in self.degree)
    is_linear = is_unbranched

    @property
    def is_unbranched(self) -> bool:
        '''Whether the monomer graph represents straight chain(s) without branching'''
        return not self.is_unbranched
    
    @property
    def terminal_monomers(self) -> Generator[int, None, None]:
        '''Generates the indices of all nodes corresponding to terminal monomers (i.e. those wiht only one outgoing bond)'''
        for node_idx, degree in self.degree:
            if degree == 1:
                yield node_idx
    termini = leaves = terminal_monomers
    

    # topological and multi-chain properties
    @property
    def num_chains(self) -> int:
        '''The number of disconnected chains represented by the MonoGraph'''
        return nx.number_connected_components(self)

    @property
    def chains(self) -> Generator['MonomerGraph', None, None]:
        '''Generates all disconnected polymers chains in the graph sequentially'''
        for cc_nodes in nx.connected_components(self):
            yield self.subgraph(cc_nodes)

    @property
    def unique_monomer_names(self) -> set[str]:
        '''The collection of unique monomer names embedded in the graph nodes'''
        return set(nx.get_node_attributes(self, self.MONOMER_NAME_ATTR).values())


    # visualization
    def draw(self, label_monomers : bool=True, label_bonds : bool=True, **kwargs) -> None: # TODO: expand arg passing (positions, matplotlib axes, etc)
        '''Visualize graph structure with NetworkX'''
        if 'pos' not in kwargs:
            kwargs['pos'] = nx.spring_layout(self) # TODO: try other layouts to see which looks best

        monomer_labels = nx.get_node_attributes(self, self.MONOMER_NAME_ATTR) if label_monomers else None
        nx.draw(self, with_labels=True, labels=monomer_labels, **kwargs)
        
        bond_labels    = nx.get_edge_attributes(self, self.BONDTYPE_ATTR) if label_bonds else None
        nx.draw_networkx_edge_labels(self, edge_labels=bond_labels, **kwargs) # TODO: add flavor labels to drawing
    visualize = draw


    # chemical information checks
    def insert_chemical_info(self, chemical_info : dict[str, dict]) -> None:
        '''Insert SMILES, SMARTS, and atom/linker count info into nodes from minimal set of monomer info templates'''
        raise NotImplemented

    def _validate(self) -> bool:
        '''Check whether the chemical information inserted into the monomer graph is valid'''
        raise NotImplemented


    # SMILES-like in-line encodings
    ## Reading string
    @classmethod
    def from_smidge_string(cls, smidge_string : str, start_node_idx : int=0) -> 'MonomerGraph':
        '''Parse a SMIDGE ("SMILES-like Monomer Interconnectivity & Degree Graph Encoding") string and read it into a networkX Graph'''
        from .smidgelib.smidgeread import SMIDGEReader

        reader = SMIDGEReader()
        return reader.read_smidge(smidge_string, start_node_idx=start_node_idx)
    from_SMIDGE = from_smidge = from_smidge_string

    ## Writing string
    def _validate_start_node_idxs(self, start_node_idxs : Optional[Union[int, Sequence[int]]]=None) -> dict[int, int]:
        '''Check if a collection of DFS traversal start indices are valid for the current graph topology'''
        # 0) if explicitly NO ids are passed, no validation is needed
        n_chains = self.num_chains
        if start_node_idxs is None:
            return {
                i : min(chain) # assign the smallest node in each component as the starting indices
                    for i, chain in enumerate(self.chains)
            }

        # 1) check that there are enough start nodes for the present number of chains
        start_node_idxs = asiterable(start_node_idxs) # convert to iterable to handle singleton values in a unified way
        n_nodes = len(start_node_idxs)
        if n_nodes != n_chains:
            quantifier = 'few' if (n_nodes < n_chains) else 'many'
            raise ValueError(f'Provided too {quantifier} chain start indices traversal of the given graph ({n_nodes} provided for {n_chains}-chain graph)')

        # 2) check that there exists a 1:1 mapping between the provided node collection and DISTINCT connected components
        cc_order_planner = DISCERNMENTSolver(nx.connected_components(self))
        if not cc_order_planner.solution_exists(start_node_idxs, unique_bins=True):
            raise ValueError('Starting node indices provided do not uniquely correspond to distinct chains')
        else:
            cc_order = next(cc_order_planner.enumerate_choices(start_node_idxs, unique_bins=True))
            return {
                chain_idx : start_node_idx
                    for (chain_idx, start_node_idx) in zip(cc_order, start_node_idxs) # the parity of this is guaranteed by the prior length match check
            }
    
    def to_smidge_string(self, start_node_idxs :  Optional[Union[int, Sequence[int]]]=None) -> str:
        '''Convert a monomer graph into a SMIDGE ("SMILES-like Monomer Interconnectivity & Degree Graph Encoding") string'''
        from .smidgelib.smidgewrite import SMIDGEWriter

        writer = SMIDGEWriter()
        return writer.write_smidge(self, start_node_idxs=start_node_idxs)
    to_smidge = to_SMIDGE = to_smidge_string

    ## Testing string translation
    def _passes_string_conversion_tests(self) -> tuple[bool, Optional[tuple[int]]]:
        '''Developer function, tests if conversion to and from graph strings preserves the graph topology invariant to the starting node
        Returns a bool of whether test passes for all possible traversal starting positions, and tuple of positions of first failure (or None if passing)'''
        for start_idxs in cartesian_product(*[chain.nodes for chain in self.chains]):
            isostr = self.to_smidge_string(start_node_idxs=start_idxs)
            isograph = self.from_smidge_string(isostr)
            if not nx.is_isomorphic(self, isograph):
                return False, start_idxs
        else:
            return True, None
MonoGraph = MonomerGraph # alias for convenience

class MonomerGraphSerializer(TypeSerializer, python_type=MonomerGraph):
    '''JSON serializer for storing MonomerGraphs as SMIDGE strings '''
    @staticmethod
    def encode(python_obj : MonomerGraph) -> str:
        return python_obj.to_smidge_string()

    @staticmethod
    def decode(json_obj : str) -> MonomerGraph:
        return MonomerGraph.from_smidge(json_obj)