'''Tools for generating and manipulating monomer connectivity graphs'''

import networkx as nx
from ..genutils.textual.delimiters import validate_braces 


class MonomerGraph(nx.Graph):
    '''A graph representation of the connectivity of monomer fragments in a polymer topology'''
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


    # chemical information checks
    def insert_chemical_info(self, monomer_info : dict[str, dict]) -> None:
        '''Insert SMILES, SMARTS, and atom/linker count info into nodes from minimal set of monomer info templates'''
        raise NotImplemented

    def _validate(self) -> bool:
        '''Check whether the chemical information inserted into the monomer graph is valid'''
        raise NotImplemented


    # SMILES-like in-line encodings
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
                monograph.add_node(curr_idx, monomer_name=mononame) #   add a new node with the current index and name
                if visited:                                 #   if previously-visited nodes exist...
                    monograph.add_edge(curr_idx, visited.pop())     #   remove the last visited node from the traversal stack and link the current node to it
                visited.append(curr_idx)                    #   add the current node to the stack of visited nodes
            elif char == '(':                               # 3) if beginning a traversal branch...
                visited.append(visited[-1])                 #   duplicate the last visited node
            elif char == ')':                               # 4) if exiting a branch...
                visited.pop()                               #   remove the last visited position and return to the previous most recent visited node
            else:                                           # 5) otherwise... 
                mononame += char                            #   must be inside a monomer block; in that case, just extend the current name
            
        return monograph

    def to_monograph_string(self : nx.Graph, start_node_idx : int=0) -> str:
        '''Convert a monomer graph into a SMILES-like monomer string'''
        neighbors = nx.dfs_successors(self, source=start_node_idx) # determine DFS ordering of nodes and their neighbors
        visited = {i : False for i in self.nodes} 
        stack = [start_node_idx]
        
        tokens = []
        while stack:
            curr_idx = stack.pop()
            if not visited[curr_idx]:                                   # 1) collect appropriate tokens for the current node, depending on whether it has already been visited
                mononame = self.nodes[curr_idx]["monomer_name"]    # get the name associated with the current monomer node (enclosed in square brackets)
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