'''Automata for reading SMIDGE strings into their graph representations'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import logging
LOGGER = logging.getLogger(__name__)

from typing import Optional
from dataclasses import dataclass, field
from abc import ABC, abstractmethod

import re

from . import MonomerGraphBondInfo, BOND_TOKEN_RE
from ..monographs import MonomerGraph
from ...genutils.textual.delimiters import validate_braces


# READER REGISTER
@dataclass
class SMIDGEReaderRegister:
    '''Information accessed by machine states during the SMIDGE writing process'''
    monograph : Optional[MonomerGraph] = field(default_factory=MonomerGraph)
    bond_info : MonomerGraphBondInfo = field(default_factory=MonomerGraphBondInfo)

    curr_token : str = field(default_factory=str)
    node_idx   : int = field(default_factory=int)
    str_buffer : str = field(default_factory=str)
    node_stack : list[int] = field(default_factory=list)


# READER STATE BASE
class SMIDGEReadState(ABC):
    '''Abstract based for reading MID graphs from SMIDGE strings'''
    @abstractmethod
    def state_action(self, register : SMIDGEReaderRegister) -> None:
        pass

    def transition(self, register : SMIDGEReaderRegister) -> 'SMIDGEReadState': # same transition rule for all states, based on the current character
        '''Define which states should follow the current one based on input'''
        READ_STATE_MAP : dict[str, SMIDGEReadState] = {
            '.' : ChainNew(),
            '[' : MonomerStart(),
            ']' : MonomerEnd(),
            '(' : BranchStart(),
            ')' : BranchEnd(),
            '<' : BondStart(),
            '>' : BondEnd(),
        }
        return READ_STATE_MAP.get(register.curr_token, Accumulate()) # treat Accumulate as default (and thereby starting) state

# CONCRETE STATES
class Accumulate(SMIDGEReadState):
    '''Collect characters into a buffer'''
    def state_action(self, register : SMIDGEReaderRegister) -> None:
        register.str_buffer += register.curr_token

class ChainNew(SMIDGEReadState):
    '''Reset actions when beginning a new chain'''
    def state_action(self, register : SMIDGEReaderRegister) -> None:
        register.node_stack.clear()

class MonomerStart(SMIDGEReadState):
    '''Begin reading in a new monomer'''
    def state_action(self, register : SMIDGEReaderRegister) -> None:
        register.str_buffer = '' # clear string buffer to read new monomer name
        register.node_idx += 1   # and increment node index

class MonomerEnd(SMIDGEReadState):
    '''Finish reading a monomer and add it to the graph'''
    def state_action(self, register : SMIDGEReaderRegister) -> None:
        register.monograph.add_node( # once a new monomer is read, add a new node corresponding to that monomer at the current index
            register.node_idx, 
            **{
                MonomerGraph.MONOMER_NAME_ATTR : register.str_buffer,
                MonomerGraph.FLAVOR_DICT_ATTR  : {}
            }
        ) #   add a new node with the current index and name
        
        if register.node_stack: # if previously-node_stack nodes exist...
            curr_node_id = register.node_idx
            prev_node_id = register.node_stack.pop() # remove the last node_stack node from the traversal stack
            register.monograph.nodes[prev_node_id][MonomerGraph.FLAVOR_DICT_ATTR][curr_node_id] = register.bond_info.incoming_flavor
            register.monograph.nodes[curr_node_id][MonomerGraph.FLAVOR_DICT_ATTR][prev_node_id] = register.bond_info.outgoing_flavor

            register.monograph.add_edge(prev_node_id, curr_node_id, **{MonomerGraph.BONDTYPE_ATTR : register.bond_info.bondtype}) # link the current node to it, with appropriate bond type
        register.node_stack.append(register.node_idx) # add the current node to the stack of node_stack nodes

class BondStart(SMIDGEReadState):
    '''Initialize reading of a bond token'''
    def state_action(self, register : SMIDGEReaderRegister) -> None:
        register.str_buffer = ''

class BondEnd(SMIDGEReadState):
    def state_action(self, register : SMIDGEReaderRegister) -> None:
        register.bond_info = MonomerGraphBondInfo.from_match(re.match(BOND_TOKEN_RE, register.str_buffer))

class BranchStart(SMIDGEReadState):
    '''Mark position of branch point in stack for backtrack'''
    def state_action(self, register : SMIDGEReaderRegister) -> None:
        register.node_stack.append(register.node_stack[-1])

class BranchEnd(SMIDGEReadState):
    '''Return to previous branch point'''
    def state_action(self, register : SMIDGEReaderRegister) -> None:
        register.node_stack.pop() # remove the last node_stack position and return to the previous most recent node_stack node
        

# READER
@dataclass
class SMIDGEReader:
    '''Pushdown automaton for translating SMIDGE strings to and from monomer graphs'''
    state : SMIDGEReadState = field(default_factory=ChainNew, init=False) # initial state is always the new chain state

    def read_smidge(self, smidge_string : str, start_node_idx : int=0) -> MonomerGraph:
        '''Parse a SMIDGE ("SMILES-like Monomer Interconnectivity & Degree Graph Encoding") string and read it into a networkX Graph'''
        validate_braces(smidge_string) # check that all braces are in order before proceeding
        register = SMIDGEReaderRegister(
            node_idx=start_node_idx-1,
        )
        for char in smidge_string:
            register.curr_token = char
            LOGGER.debug(f'Current {register.__class__.__name__}: {register!r}')
            LOGGER.debug(f'Current {self.__class__.__name__} state: {self.state.__class__.__name__}')
            self.state = self.state.transition(register)
            self.state.state_action(register) # TODO: figure out how to make this work with the order of actions and transitions reversed

        return register.monograph