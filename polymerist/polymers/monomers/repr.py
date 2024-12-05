'''For representing monomer information'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Generator, Optional, TypeAlias, Union
from dataclasses import dataclass, field

from collections import defaultdict
from rdkit import Chem
from rdkit.Chem.rdchem import Mol

from ...genutils.iteration import iter_len
from ...genutils.fileutils.jsonio.jsonify import make_jsonifiable
from ...rdutils.bonding.portlib import get_num_ports


ResidueSmarts : TypeAlias = dict[str, list[str]] # monomer SMARTS strings keyed by residue name

# MAIN REPRESENTATION CLASS
@make_jsonifiable
@dataclass
class MonomerGroup:
    '''Stores collections of residue-labelled monomer SMARTS'''
    monomers : ResidueSmarts = field(default_factory=dict)
    term_orient : dict[str, str] = field(default_factory=dict) # keys are either "head" or "tail", values are the names of residues in "monomers"

    @staticmethod
    def is_terminal(monomer : Mol) -> bool:
        '''Determine whether or not a monomer is terminal'''
        return get_num_ports(monomer) == 1

    # ATTRIBUTE PROPERTIES AND ALIASES
    @property
    def SMARTS(self) -> ResidueSmarts:
        '''Alias of legacy "monomers" attribute'''
        return self.monomers # alias of legacy name for convenience
    
    def iter_rdmols(self, term_only : Optional[bool]=None) -> Generator[tuple[str, Mol], None, None]:
        '''
        Generate (residue name, RDKit Mol) pairs of all monomers present
        Simplifies iteration over internal lists of monomer Mols

        Can optionally filter by monomer termination:
            term_only=True  -> only terminal monomers
            term_only=False -> only middle monomers
            term_only=None  -> all monomers
        '''
        for resname, SMARTS_list in self.monomers.items():
            for SMARTS in SMARTS_list:
                monomer = Chem.MolFromSmarts(SMARTS)
                if (term_only is None) or (MonomerGroup.is_terminal(monomer) == term_only):
                    yield (resname, monomer)

    def rdmols(self, term_only : Optional[bool]=None) -> dict[str, list[Mol]]:
        '''
        Returns dict of RDKit Mol lists keyed by residue name

        Can optionally filter by monomer termination:
            term_only=True  -> only terminal monomers
            term_only=False -> only middle monomers
            term_only=None  -> all monomers
        '''
        rdmol_dict = defaultdict(list)
        for resname, rdmol in self.iter_rdmols(term_only=term_only):
            rdmol_dict[resname].append(rdmol)

        return rdmol_dict
    
    @property
    def n_monomers(self) -> int:
        '''Returns number of present monomers
        Multiple monomers with the same residue name are considered distinct'''
        return iter_len(self.iter_rdmols(term_only=None))
    
    # VALIDATION AND PROPERTY CHECKS
    @property
    def _is_valid(self) -> bool:
        '''Check that types and formatting are correct'''
        for resname, SMARTS_list in self.monomers.items():
            if not (isinstance(resname, str) and isinstance(SMARTS_list, list)):
                return False
        else:
            return True # valid only if none of the SMARTS lists fail
        
    @property
    def has_valid_linear_term_orient(self) -> bool:
        '''Check whether terminal group orientations are sufficient to define a linear polymer'''
        return (
            bool(self.term_orient)                                         # check that: 1) term group orientations are non-empty...
            and set(self.term_orient.keys()) == {'head', 'tail'}                       # 2) ...orientation labels are only "head" and "tail" (in any order)...
            and all(resname in self.monomers for resname in self.term_orient.values()) # 3) ... and all term group keys match a present monomer
        )
    
    # COMPOSITION AND I/O METHODS
    def __add__(self, other : 'MonomerGroup') -> 'MonomerGroup':
        '''Content-aware method of merging multiple sets of monomer info via the addition operator'''
        cls = self.__class__
        if not isinstance(other, cls):
            raise NotImplementedError(f'Can only merge {cls.__name__} with another {cls.__name__}, not object of type {type(other)}')

        return MonomerGroup(monomers={**self.monomers, **other.monomers})

    __radd__ = __add__ # support reverse addition

    # CHEMICAL INFORMATION
    def unique(self, cap_group : Union[str, Mol]=Chem.MolFromSmarts('[H]-[*]')) -> 'MonomerGroup':
        '''Return a MonomerGroup containing only the unique monomers present, given a particular port saturating group (by default just a hydrogen)'''
        raise NotImplementedError
        # unique_mono = set()
        # for SMARTS in monomer_smarts.values():
        #     monomer = Chem.MolFromSmarts(SMARTS)
        #     clear_atom_map_nums(monomer, in_place=True) 
        #     hydrogenate_monomer_ports(monomer, in_place=True) 
        #     unique_mono.add(Chem.MolToSmiles(monomer)) # TODO : eventually make this SMARTS-based (or even better RDKit Mol-based); can't for now since hydrogenated fragments don't equate

        # return unique_mono

    def is_homopolymer(self) -> bool:
        '''Identify if a polymer is a homopolymer (i.e. only 1 type of middle monomer)'''
        # n_mid, n_term = count_middle_and_term_mono(monomer_smarts) # TODO : reimplement with comparison of port-hydrogenated monomers
        # return (n_mid == 1)
        return (len(self.unique()) == 1) # by definition, a homopolymer only has 1 unique class of monomer

    # GRAPH INFORMATION
    @property
    def is_branchable(self) -> bool:
        '''Whether it is possible to generate a branched polymer from this set of monomers'''
        return any(
            get_num_ports(monomer) > 2
                for (resname, monomer) in self.iter_rdmols(term_only=None)
        )
    
    @property
    def is_linear(self) -> bool:
        '''Whether a group of monomers can ONLY be assembled into a linear chain'''
        return not self.is_branchable

    @property
    def is_linear_homopolymer(self) -> bool:
        '''Identify if a polymer is a linear homopolymer'''
        return self.is_linear and self.is_homopolymer

    @property
    def num_mid_and_term(self) -> tuple[int, int]:
        '''Counts of how many of the monomers are middle vs terminal, respectively'''
        group_counts = [0, 0]
        for (resname, monomer) in self.iter_rdmols(term_only=None): # TODO : consider reimplementing using new term group filtering option
            group_counts[self.is_terminal(monomer)] += 1 # index by bool
        
        return tuple(group_counts) # convert to tuple