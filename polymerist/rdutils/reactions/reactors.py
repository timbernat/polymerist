'''Classes for implementing reactions with respect to some set of reactant RDKit Mols'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Generator, Sequence
from dataclasses import dataclass, field

from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdmolops import SanitizeMol, SanitizeFlags, SANITIZE_ALL

from .reactexc import ReactantTemplateMismatch
from .reactions import AnnotatedReaction
from .fragment import IBIS, ReseparateRGroups

from ..rdprops.atomprops import clear_atom_props
from ..rdprops.bondprops import clear_bond_props
from ..chemlabel import clear_atom_map_nums


@dataclass
class PolymerizationReactor:
    '''Reactor which exhaustively generates monomers fragments according to a given a polymerization mechanism'''
    rxn_schema : AnnotatedReaction
    fragment_strategy : IBIS = field(default_factory=ReseparateRGroups)
    
    def propagate(
        self,
        monomers : Sequence[Mol],
        clear_map_nums : bool=True,
        sanitize_ops : SanitizeFlags=SANITIZE_ALL,
     ) -> Generator[tuple[list[Mol], list[Mol]], None, None]:
        '''Keep reacting and fragmenting a pair of monomers until all reactive sites have been reacted
        Returns fragment pairs at each step of the chain propagation process'''
        reactants = monomers # initialize reactive pair with monomers
        while True: # check if the reactants can be applied under the reaction template
            try:
                adducts = self.rxn_schema.react(reactants, repetitions=1, sanitize_ops=sanitize_ops) # can't clear properties yet, otherwise intermonomer bond finder would have nothing to work with
            except ReactantTemplateMismatch:
                break
            
            fragments : list[Mol] = []
            for product in adducts: # DEVNOTE: consider doing fragmentation on the combined molecule made up of all products?
                for fragment in self.fragment_strategy.produce_fragments(product, separate=True):
                    clear_atom_props(fragment, in_place=True) # essential to avoid reaction mapping info from prior steps from contaminating future ones
                    clear_bond_props(fragment, in_place=True)
                    SanitizeMol(fragment, sanitizeOps=sanitize_ops)
                    fragments.append(fragment)

                if clear_map_nums: # NOTE : CRITICAL that this be done after fragmentation step, which RELIES on map numbers being present
                    clear_atom_map_nums(product, in_place=True)
                    
            yield adducts, fragments # yield the adduct Mol and any subsequent resulting reactive fragments
            reactants = fragments # set fragments from current round of polymerization as reactants for next round
            
    def propagate_iterative(self) -> None:
        ...