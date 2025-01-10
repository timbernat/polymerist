'''For representing monomer information'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import logging
LOGGER = logging.getLogger(__name__)

from typing import Generator, Optional, Iterable, Union
from dataclasses import dataclass, field

from itertools import cycle
from collections import defaultdict

from rdkit import Chem

from ...genutils.iteration import iter_len
from ...genutils.fileutils.jsonio.jsonify import make_jsonifiable

from ...smileslib.primitives import Smarts, is_valid_SMARTS
from ...rdutils.bonding.portlib import get_num_ports


# MAIN REPRESENTATION CLASS
@make_jsonifiable
@dataclass
class MonomerGroup:
    '''Stores collections of residue-labelled monomer SMARTS'''
    monomers : dict[str, Union[Smarts, list[Smarts]]] = field(default_factory=dict)
    term_orient : dict[str, str] = field(default_factory=dict) # keys are either "head" or "tail", values are the names of residues in "monomers"

    # MONOMER ADDITION AND VALIDATION
    def __post_init__(self) -> None:
        # Encase bare SMARTS into lists and check that all monomer SMARTS are valid
        monomers_init = self.monomers # store inputted values
        self.monomers = {} # clear monomers and re-add one-at-a-time
        for resname, smarts in monomers_init.items():
            self.add_monomer(resname, smarts)
        # DEV: opted to forgo term_orient check for now, as modifying this violates the read-only data model aimed for here
                
    def _add_monomer(self, resname : str, smarts : Smarts) -> None:
        '''Add a new monomer to the templates already stored within, subject to validation checks'''
        if not isinstance(smarts, str): 
            raise TypeError(f'Values of monomers must be either SMARTS strings or lists of SMARTS strings, not "{type(smarts).__name__}"')
        # DEV: include check for empty string? (technically still a valid SMARTS string, but a pretty pathological one at that)
        if not is_valid_SMARTS(smarts):
            raise ValueError(f'Provided invalid monomer SMARTS string for {resname}: "{smarts}"') 
        # DEV: decide whether or not SMILES expansion and spec-compliance should be enforced here or shunted off to the user 
        
        if resname in self.monomers:
            existing_resgroup = self.monomers[resname]
            if isinstance(existing_resgroup, list) and (smarts not in existing_resgroup):
                LOGGER.debug(f'Extending existing residue category "{resname}" with SMARTS {smarts}')
                self.monomers[resname].append(smarts)
        else:
            LOGGER.debug(f'Creating new residue category "{resname}", containing singular SMARTS ["{smarts}"])')
            self.monomers[resname] = [smarts]
            
    def _add_monomers(self, resname : str, smarts_container : Iterable[Smarts]) -> None:
        '''Add new monomers to the templates already stored within, subject to validation checks, from an iterable container'''
        for smarts in smarts_container:
            self._add_monomer(resname, smarts)
    
    def add_monomer(self, resname : str, smarts : Union[Smarts, Iterable[Smarts]]) -> None:
        '''Register new monomers, either directly from SMARTS or from a container of SMARTS'''
        if isinstance(smarts, Iterable) and not isinstance(smarts, str): # don;t want to insert one character at a time if a string is in fact provided
            self._add_monomers(resname, smarts)
        else:
            self._add_monomer(resname, smarts) # assume any other inputs are singular values or strings 
    
    # DUNDER "MAGIC" METHODS
    def __getitem__(self, resname : str) -> str:
        '''Convenience method to access .monomers directly from instance'''
        return self.monomers[resname] # NOTE: deliberately avoid "get()" here to propagate KeyError
        # BUG: user can directly append to the returned value to forgo monomer validation checks;
        # this is not unit to __getitem__ but rather a consequence of thinly-wrapping builtin types

    def __setitem__(self, resname : str, smarts : Smarts) -> str:
        '''Convenience method to access .monomers directly from instance'''
        self.add_monomer(resname, smarts)
        
    def __hash__(self) -> int:
        '''Hash based on monomer SMARTS and terminal orientation in a canonical order'''
        # TOSELF: this is far from bulletproof, viz. canonicalzation of SMARTS, list value sorting, etc
        return hash(f'{sorted(self.monomers.items())}{sorted(self.term_orient.items())}')
    
    # ATTRIBUTE PROPERTIES AND ALIASES
    @staticmethod
    def is_terminal(monomer : Chem.Mol) -> bool:
        '''Determine whether or not a monomer is terminal'''
        return get_num_ports(monomer) == 1
    
    @property
    def SMARTS(self) -> dict[str, list[Smarts]]:
        '''Alias of legacy "monomers" attribute'''
        return self.monomers # alias of legacy name for convenience
    
    # ITERATION OVER STORED MOLECULE FRAGMENTS
    def iter_rdmols(self, term_only : Optional[bool]=None) -> Generator[tuple[str, Chem.Mol], None, None]:
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

    def rdmols(self, term_only : Optional[bool]=None) -> dict[str, list[Chem.Mol]]:
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
    
    def contributions(self, term_only : Optional[bool]=None) -> dict[str, list[int]]:
        '''Returns dict of the number of real (i.e. non-linker) atoms in each residue list'''
        return {
            resname : [mol.GetNumAtoms() - get_num_ports(mol) for mol in mol_list]
                for resname, mol_list in self.rdmols(term_only=term_only).items()
        }
    
    @property
    def n_monomers(self) -> int:
        '''Returns number of present monomers; multiple monomers under the same residue name are considered distinct'''
        return iter_len(self.iter_rdmols(term_only=None))
    
    # END GROUP DETERMINATION 
    def linear_end_groups(self) -> dict[str, tuple[str, Chem.Mol]]:
        '''
        Returns head-and-tail end group residue names and Mol objects as defined by term_orient
        
        If term orient is undefined, will automatically take then first 
        <= 2 terminal groups available to be the end groups
        
        Returns
        -------
        end_groups : dict[str, tuple[str, Chem.Mol]]
            A dict whose keys are any of {'head', 'tail'} and whose
            values are 2-tuples of residue names and Mols for the corresponding monomer
        '''
        if self.term_orient and set(self.term_orient.keys()) == {'head', 'tail'}:
            LOGGER.info(f'Using user-defined terminal group orientation {self.term_orient}')
            monomer_iters = {
                resname : cycle(smarts_list) 
                    for resname, smarts_list in self.rdmols(term_only=True).items()
            } # cycle handles degenerate end group case correctly
            
            return {
                head_or_tail : (resname, next(monomer_iters[resname])) # will raise KeyError if any of the resnames are not present
                    for head_or_tail, resname in self.term_orient.items()
            }
        else:
            term_orient_auto : dict[str, Smarts] = {}
            end_groups_auto  : dict[str, Chem.Mol] = {}
            for head_or_tail, (resname, rdmol) in zip(['head', 'tail'], self.iter_rdmols(term_only=True)): # zip will bottom out early if fewer than 2 terminal monomers are present
                term_orient_auto[head_or_tail] = resname # populate purely for logging
                end_groups_auto[head_or_tail]  = (resname, rdmol)
            LOGGER.warning(f'No valid terminal monomer orientations defined, auto-assigned orientations "{term_orient_auto}"; USER SHOULD VERIFY THIS YIELDS A CHEMICALLY-VALID POLYMER!')
                
            return end_groups_auto
    
    # COMPOSITION METHODS
    def __add__(self, other : 'MonomerGroup') -> 'MonomerGroup':
        '''Content-aware method of merging multiple sets of monomer info via the addition operator'''
        cls = self.__class__
        if not isinstance(other, cls):
            raise NotImplementedError(f'Can only merge {cls.__name__} with another {cls.__name__}, not object of type {type(other)}')
        # TODO: figure out how to handle combination of term group orientation gracefully (ignoring for now)
        return MonomerGroup(monomers={**self.monomers, **other.monomers})

    __radd__ = __add__ # support reverse addition

    # CHEMICAL INFORMATION
    def unique(self, cap_group : Union[Smarts, Chem.Mol]=Chem.MolFromSmarts('[H]-[*]')) -> 'MonomerGroup':
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