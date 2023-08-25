'''For representing monomer information'''

from typing import Generator, TypeAlias, Union

from dataclasses import dataclass, field
from rdkit import Chem

from ..genutils.fileutils.jsonio import JSONifiable
from ..rdutils.rdtypes import RDMol
from ..rdutils.amalgamation.portlib import get_num_ports


ResidueSmarts : TypeAlias = dict[str, list[str]] # monomer SMARTS strings keyed by residue name


# MAIN REPRESENTATION CLASS
@dataclass
class MonomerGroup(JSONifiable):
    '''Stores collections of residue-labelled monomer SMARTS'''
    monomers : ResidueSmarts = field(default_factory=dict)

    @staticmethod
    def is_terminal(monomer : RDMol) -> bool:
        '''Determine whether or not a monomer is terminal'''
        return get_num_ports(monomer) == 1

    # ATTRIBUTE PROPERTIES AND ALIASES
    @property
    def SMARTS(self) -> ResidueSmarts:
        '''Alias of legacy "monomers" attribute'''
        return self.monomers # alias of legacy name for convenience
    
    @property
    def rdmols(self) -> dict[str, list[RDMol]]:
        '''Dict of RDKit Mols produced by ResidueSmarts, for convenience'''
        return {
            resname : [Chem.MolFromSmarts(SMARTS) for SMARTS in SMARTS_list]
                for resname, SMARTS_list in self.monomers.items()
        }

    @property
    def iter_rdmols(self) -> Generator[tuple[str, RDMol], None, None]:
        '''Simplifies iteration over internal lists'''
        for resname, monomer_list in self.rdmols.items():
            for monomer in monomer_list:
                yield (resname, monomer)

    @property
    def _is_valid(self) -> bool:
        '''Check that types and formatting are correct'''
        for resname, SMARTS_list in self.monomers.items():
            if not (isinstance(resname, str) and isinstance(SMARTS_list, list)):
                return False
        else:
            return True
    
    # COMPOSITION AND I/O METHODS
    def __add__(self, other : 'MonomerGroup') -> 'MonomerGroup':
        '''Content-aware method of merging multiple sets of monomer info via the addition operator'''
        cls = self.__class__
        if not isinstance(other, cls):
            raise NotImplementedError(f'Can only merge {cls.__name__} with another {cls.__name__}, not object of type {type(other)}')

        return MonomerGroup(monomers={**self.monomers, **other.monomers})

    __radd__ = __add__ # support reverse addition

    # CHEMICAL INFORMATION
    def unique(self, cap_group : Union[str, RDMol]=Chem.MolFromSmarts('[H]-[*]')) -> 'MonomerGroup':
        '''Return a MonomerGroup containing only the unique monomers present, given a particular port saturating group (by default just a hydrogen)'''
        raise NotImplemented
        # unique_mono = set()
        # for SMARTS in monomer_smarts.values():
        #     monomer = Chem.MolFromSmarts(SMARTS)
        #     clear_atom_map_nums(monomer, in_place=True) 
        #     hydrogenate_monomer_ports(monomer, in_place=True) 
        #     unique_mono.add(Chem.MolToSmiles(monomer)) # TODO : eventually make this SMART-based (or even better RDMol-based); can't for now since hydrogenated fragments don't equate

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
                for (resname, monomer) in self.iter_rdmols
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
        for (resname, monomer) in self.iter_rdmols:
            group_counts[self.is_terminal(monomer)] += 1 # index by bool
        
        return tuple(group_counts) # convert to tuple