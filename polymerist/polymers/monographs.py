'''Tools for generating and manipulating monomer connectivity graphs'''

from typing import ClassVar, Generator, Optional, Sequence, Union

import networkx as nx
from itertools import product as cartesian_product

from ..genutils.iteration import asiterable
from ..genutils.sequences import bin_ids_forming_sequence
from ..genutils.textual.delimiters import validate_braces 


class MonomerGraph(nx.Graph):
    '''A graph representation of the connectivity of monomer fragments in a polymer topology'''
    MONOMER_NAME_ATTR : ClassVar[str] = 'monomer_name' # node attribute name to assign the name of each monomer to

    # connectivity properties
    @property
    def num_monomers(self) -> int:
        '''Number of monomer units represented in the current polymer'''
        return self.number_of_nodes()
    DOP = num_monomers

    @property
    def is_unbranched(self) -> bool:
        '''Whether the monomer graph represents straight chain(s) without branching'''
        return all(node_deg <= 2 for node_id, node_deg in self.degree)
    is_linear = is_unbranched

    @property
    def is_unbranched(self) -> bool:
        '''Whether the monomer graph represents straight chain(s) without branching'''
        return not self.is_unbranched
    

    # topological and multi-chain settings
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
        '''The collection of unique monomer '''
        return set(nx.get_node_attributes(self, self.MONOMER_NAME_ATTR).values())


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
    def from_monograph_string(cls, monostring : str, start_node_idx : int=0) -> 'MonomerGraph':
        '''Parse a SMILES-like monomer graph string and read it into a networkX Graph'''
        validate_braces(monostring) # check that all braces are in order before proceeding
        visited : list[int] = []
        curr_idx = start_node_idx - 1
        mononame = ''
        
        monograph = cls() # create empty instance with class initializer
        for char in monostring:
            if char == '[':                                 # 1) if reached a new monomer block...
                mononame = ''                               #   clear the current monomer name
                curr_idx += 1                               #   and advance the current position index
            elif char == ']':                               # 2) if reached the end of a monomer block...
                monograph.add_node(curr_idx, **{cls.MONOMER_NAME_ATTR : mononame}) #   add a new node with the current index and name
                if visited:                                 #   if previously-visited nodes exist...
                    monograph.add_edge(curr_idx, visited.pop())     #   remove the last visited node from the traversal stack and link the current node to it
                visited.append(curr_idx)                    #   add the current node to the stack of visited nodes
            elif char == '(':                               # 3) if beginning a traversal branch...
                visited.append(visited[-1])                 #   duplicate the last visited node
            elif char == ')':                               # 4) if exiting a branch...
                visited.pop()                               #   remove the last visited position and return to the previous most recent visited node
            elif char == '.':                               # 5) if reaching the end of a connected component... | TOSELF: opted for this over str.split('.') as this already tracks node sequence
                visited.clear()                             #   clear the visited stack to restart the traversal for the new component (ensures it is not connected to the previous one)
            else:                                           # 6) otherwise... 
                mononame += char                            #   must be inside a monomer block; in that case, just extend the current name
            
        return monograph

    ## Writing string
    def _validate_start_node_idxs(self : nx.Graph, start_node_idxs : Optional[Union[int, Sequence[int]]]=None) -> dict[int, int]:
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
        possible_cc_orders = bin_ids_forming_sequence(start_node_idxs, nx.connected_components(self), unique_bins=True)
        try:
            cc_order = next(possible_cc_orders)
            return {
                chain_idx : start_node_idx
                    for (chain_idx, start_node_idx) in zip(cc_order, start_node_idxs) # the parity of this is guaranteed by the prior length match check
        }
        except StopIteration:
            raise ValueError('Starting node indices provided do not uniquely correspond to distinct chains')

    def _chain_to_monograph_string(self : nx.Graph, start_node_idx : Optional[int]=None) -> str:
        '''Convert an individual chain in monomer graph into a SMILES-like monomer string'''
        neighbors = nx.dfs_successors(self, source=start_node_idx) # determine DFS ordering of nodes and their neighbors
        visited = {i : False for i in self.nodes} 
        stack = [start_node_idx]
        
        tokens = []
        while stack:
            curr_idx = stack.pop()
            if not visited[curr_idx]:                                   # 1) collect appropriate tokens for the current node, depending on whether it has already been visited
                mononame = self.nodes[curr_idx][self.MONOMER_NAME_ATTR] # get the name associated with the current monomer node (enclosed in square brackets)
                tokens.append(f'[{mononame}]')                          #   push the current node's monomer label onto the result stack
                visited[curr_idx] = True                                #   and mark as having been visited
            else:                                                       # otherwise, if returning to an already-visited node...
                tokens.append(')')                                      #   close the branch that must have led to this node
                
            if (successors := neighbors.get(curr_idx, [])):             # 2) get the remaining DFS successors of the current node, proceeding with checks if nonempty...
                if (len(successors) > 1):                               # if multiple unvisited successors are present...
                    tokens.append('(')                                  #   initialize a new branch point indicator 
                    stack.append(curr_idx)                              #   and mark the current node 
                stack.append(successors.pop(0))                         # push the first available successor node to be visited next

        return ''.join(tokens)
    
    def to_monograph_string(self : nx.Graph, start_node_idxs :  Optional[Union[int, Sequence[int]]]=None) -> str:
        '''Convert a monomer graph into a SMILES-like monomer string'''
        chain_starts = self._validate_start_node_idxs(start_node_idxs)
        return '.'.join(
            self._chain_to_monograph_string(start_node_idx=chain_starts[i])
                for i in range(self.num_chains)
        )

    ## Testing string translation
    def _passes_string_conversion_tests(self) -> tuple[bool, Optional[tuple[int]]]:
        '''Developer function, tests if conversion to and from graph strings preserves the graph topology invariant to the starting node
        Returns a bool of whether test passes for all possible traversal starting positions, and tuple of positions of first failure (or None if passing)'''
        for start_idxs in cartesian_product(*[chain.nodes for chain in self.chains]):
            isostr = self.to_monograph_string(start_node_idxs=start_idxs)
            isograph = self.from_monograph_string(isostr)
            if not nx.is_isomorphic(self, isograph):
                return False, start_idxs
        else:
            return True, None