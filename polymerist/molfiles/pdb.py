'''PDB file formatting tools'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from dataclasses import dataclass, field
from collections import Counter


@dataclass(frozen=True)
class SerialAtomLabeller:
    '''
    For assigning unique numbered atom names based on their
    order of appearance within a molecule and elemental class
    
    Useful, for example, in generating unique atom names for a PDB file
    
    Parameters
    ----------
    atom_label_width : int , default 4      
        Exact length alloted for any generated atom label
        Labels shorter than this are right-padded with spaces,
        while labels longer than this are truncated
        
        Default of 4 is the chosen to be compatible with the PDB specification ("Atom name: lines 13-16, left-justified")
        https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
    include_elem_idx : bool, default True  
        Whether to attach a numerical element-index postfix to atom labels
        
        E.g. with atom_label_width=4, the fifth carbon in a topology  
        will be labelled as "C004" with include_elem_idx=True, 
        while labelled as "C   " with include_elem_idx=False, 
    default_elem_idx : int, default 0
        Starting index for each element category
        By default, is 0-indexed; MUST BE POSITIVE
    '''
    atom_label_width : int = 4
    include_elem_idx : bool = True
    default_elem_idx : int = 0
    
    element_counter : Counter = field(init=False, default_factory=Counter)
    
    def __post_init__(self) -> None:
        '''Check ranges on input values'''
        if self.atom_label_width < 0:
            raise ValueError(f'Must provide a non-negative number of index digits to include (provided {self.atom_label_width})')

        if self.default_elem_idx < 0:
            raise ValueError(f'Must provide a non-negative starting index for element indices (provided {self.default_elem_idx})')
    
    def get_atom_label(self, elem_symbol : str) -> str:
        '''
        Obtain a numbered atom label for an atom based on its element, 
        updating the underlying element context in the process
        '''
        if not isinstance(elem_symbol, str):
            raise TypeError(f'Must pass symbol of atom\'s element as str (not type {type(elem_symbol).__name__})')
        
        if elem_symbol not in self.element_counter: # initialize first occurence to starting value
            self.element_counter[elem_symbol] = self.default_elem_idx
            
        atom_idx_label : str = ''
        if self.include_elem_idx:
            atom_idx = self.element_counter[elem_symbol]
            num_idx_digits = max(self.atom_label_width - len(elem_symbol), 0) # number of symbols left over for an atom index
            atom_idx_label = f'{atom_idx:0{num_idx_digits}d}'
        
        atom_name = f'{elem_symbol}{atom_idx_label}'
        atom_name = atom_name.ljust(self.atom_label_width, ' ')[:self.atom_label_width] # pad with spaces if too short, or truncate if too long
        assert(len(atom_name) <= self.atom_label_width) # perfunctory check to make sure things are working as expected
        
        self.element_counter[elem_symbol] += 1 # update tally with addition of new occurence of a particular element
        
        return atom_name
    