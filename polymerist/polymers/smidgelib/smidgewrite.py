'''Automata for reading SMIDGE strings into their graph representations'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import logging
LOGGER = logging.getLogger(__name__)

from typing import Optional, Sequence, Union
from dataclasses import dataclass, field
from abc import ABC, abstractmethod

import networkx as nx

from .smidgebonds import MonomerGraphBondInfo
from ..monographs import MonomerGraph


@dataclass
class SMIDGEWriterRegister:
    '''Information accessed by machine states during the SMIDGE writing process'''
    monograph    : MonomerGraph         = field(default_factory=MonomerGraph) # this is the only initializable field
    bond_info    : MonomerGraphBondInfo = field(default_factory=MonomerGraphBondInfo)
    tokens       : list[int]            = field(default_factory=list)

    curr_idx     : int                  = field(default_factory=int)
    visited      : set                  = field(default_factory=set)
    node_stack   : list[int]            = field(default_factory=list)
    successors   : dict[int, list[int]] = field(default_factory=dict)
    predecessors : dict[int, int]       = field(default_factory=dict)

    finished     : bool                 = field(default_factory=bool)

    @property
    def prev_idx(self) -> Optional[int]:
        '''The index of the precursor to the currently active node'''
        return self.predecessors.get(self.curr_idx)
    
    @property
    def node_successors(self) -> list[int]:
        '''The nodes which follow the current node in the DFS traversal (could be empty)'''
        return self.successors.get(self.curr_idx)


# STATE BASE
class SMIDGEWriterState(ABC):
    '''Abstract base for automaton states used in writing MID graphs to SMIDGE strings'''
    # NOTE: state_action() is deliberately NOT an abstract method, as doing nothing while in a state is perfectly valid
    def state_action(self, register : SMIDGEWriterRegister) -> 'SMIDGEWriterState':
        pass

    @abstractmethod
    def transition(self, register : SMIDGEWriterRegister) -> 'SMIDGEWriterState':
        '''Define which states should follow the current one based on input'''
        pass

# CONCRETE STATES
class CheckIsVisitedNode(SMIDGEWriterState):
    '''Read a node index off the stack and check if it has already been visited'''
    def state_action(self, register: SMIDGEWriterRegister) -> SMIDGEWriterState:
        register.curr_idx = register.node_stack.pop()

    def transition(self, register: SMIDGEWriterRegister) -> SMIDGEWriterState:
        if register.curr_idx in register.visited:
            return DetermineNodeSuccessors()
        else:
            return VisitNode()

class VisitNode(SMIDGEWriterState):
    '''Visit a node, check if it has predecessors, and mark it as visited'''
    def state_action(self, register: SMIDGEWriterRegister) -> SMIDGEWriterState:
        if register.prev_idx is not None:
            register.bond_info = MonomerGraphBondInfo(
                incoming_flavor=register.monograph.nodes[register.prev_idx].get(MonomerGraph.FLAVOR_DICT_ATTR, {}).get(register.curr_idx),
                bondtype=register.monograph.edges[register.prev_idx, register.curr_idx].get(MonomerGraph.BONDTYPE_ATTR),
                outgoing_flavor=register.monograph.nodes[register.curr_idx].get(MonomerGraph.FLAVOR_DICT_ATTR, {}).get(register.prev_idx), # NOTE : the same as incoming_flavor, but with the order of nodes reversed
            )
        else:
            register.bond_info = MonomerGraphBondInfo()
        register.visited.add(register.curr_idx) # mark the current node as visited

    def transition(self, register: SMIDGEWriterRegister) -> SMIDGEWriterState:
        if register.prev_idx is not None:
            return WriteEdge()
        else:
            return WriteNode()

class WriteNode(SMIDGEWriterState):
    '''Add an entry for the current node to the output tokens'''
    def state_action(self, register: SMIDGEWriterRegister) -> SMIDGEWriterState:
        monomer_token = f'[{register.monograph.nodes[register.curr_idx][register.monograph.MONOMER_NAME_ATTR]}]'
        register.tokens.append(monomer_token)
        
        if not register.successors:
            register.finished = True # will implicitly cause this to be the final machine state when parsin
    
    def transition(self, register: SMIDGEWriterRegister) -> SMIDGEWriterState:
        return DetermineNodeSuccessors()

class WriteEdge(SMIDGEWriterState):
    '''Add an entry for the current edge to the output tokens'''
    def state_action(self, register: SMIDGEWriterRegister) -> SMIDGEWriterState:
        register.tokens.append(f'<{register.bond_info!s}>')

    def transition(self, register: SMIDGEWriterRegister) -> SMIDGEWriterState:
        return WriteNode()

class DetermineNodeSuccessors(SMIDGEWriterState):
    '''Determine branching and next-node behavior'''
    def transition(self, register: SMIDGEWriterRegister) -> SMIDGEWriterState:
        if register.node_successors:
            if len(register.node_successors) > 1:
                return BranchStart()
            else:
                return PushNodeSuccessor()
        else:
            if register.successors:
                return BranchEnd() # always return a branch end condition
            else:
                return UpdateNodeSuccessors() # EXCEPT when the current node is also the last node

class PushNodeSuccessor(SMIDGEWriterState):
    '''Add a nodes DFS successor to the traversal stack'''
    def state_action(self, register: SMIDGEWriterRegister) -> SMIDGEWriterState:
        register.node_stack.append(register.node_successors.pop(0))
    
    def transition(self, register: SMIDGEWriterRegister) -> SMIDGEWriterState:
        return UpdateNodeSuccessors()

class UpdateNodeSuccessors(SMIDGEWriterState):
    '''Remove any node with no successors from the DFS traversal map'''
    def state_action(self, register: SMIDGEWriterRegister) -> SMIDGEWriterState:
        if not register.node_successors:
            register.successors.pop(register.curr_idx, None) # remove the current node from the map of possible successors (None default silences KeyError if only implicitly present)

    def transition(self, register: SMIDGEWriterRegister) -> SMIDGEWriterState:
        return CheckIsVisitedNode()

class BranchStart(SMIDGEWriterState):
    '''Mark position of branch point in stack for backtrack'''
    def state_action(self, register: SMIDGEWriterRegister) -> SMIDGEWriterState:
        register.node_stack.append(register.curr_idx) # push the current node onto the stack again to mark for revisiting later
        register.tokens.append('(')

    def transition(self, register: SMIDGEWriterRegister) -> SMIDGEWriterState:
        return PushNodeSuccessor()

class BranchEnd(SMIDGEWriterState):
    '''Mark the end of a branch'''
    def state_action(self, register: SMIDGEWriterRegister) -> SMIDGEWriterState:
        register.tokens.append(')')

    def transition(self, register: SMIDGEWriterRegister) -> SMIDGEWriterState:
        return UpdateNodeSuccessors()


# WRITER
@dataclass
class SMIDGEWriter:
    '''Register machine for translating SMIDGE strings to and from monomer graphs'''
    state : SMIDGEWriterState = field(default_factory=CheckIsVisitedNode, init=False) # initial state is always the new chain state

    def write_smidge_connected_component(self, monograph : MonomerGraph, start_node_idx : Optional[int]=None) -> MonomerGraph:
        '''Parse a single connected component of a monomer graph into a SMIDGE string'''
        self.state = CheckIsVisitedNode()
        register = SMIDGEWriterRegister(
            monograph=monograph,
            curr_idx=start_node_idx,
            successors=nx.dfs_successors(monograph, source=start_node_idx),
            predecessors=nx.dfs_predecessors(monograph, source=start_node_idx),
            node_stack=[start_node_idx],
            finished=False
        )

        while not register.finished:
            LOGGER.debug(f'Current {register.__class__.__name__}: {register!r}')
            LOGGER.debug(f'Current {self.__class__.__name__} state: {self.state.__class__.__name__}')
            self.state.state_action(register)
            self.state = self.state.transition(register)
        return ''.join(register.tokens)
    
    def write_smidge(self, monograph : MonomerGraph, start_node_idxs :  Optional[Union[int, Sequence[int]]]=None) -> str:
        '''Parse all connected components of a monograph into a complete SMIDGE string'''
        chain_starts = monograph._validate_start_node_idxs(start_node_idxs)
        return '.'.join(
            self.write_smidge_connected_component(chain, start_node_idx=chain_starts[i])
                for i, chain in enumerate(monograph.chains)
        )