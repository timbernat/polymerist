'''Classes for implementing reactions with respect to some set of reactant RDKit Mols'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import logging
LOGGER = logging.getLogger(__name__)

from typing import Generator, Sequence
from dataclasses import dataclass, field

from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdmolops import SanitizeMol, SanitizeFlags, SANITIZE_ALL
from rdkit.Chem.rdmolops import AromaticityModel, AROMATICITY_RDKIT, AROMATICITY_MDL

from .reactions import AnnotatedReaction
from .fragment import IBIS, ReseparateRGroups

from ..rdprops.atomprops import clear_atom_props
from ..rdprops.bondprops import clear_bond_props
from ..chemlabel import clear_atom_map_nums
from ..sanitization import sanitize

from ...smileslib.cleanup import canonical_SMILES_from_mol


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
        aromaticity_model : AromaticityModel=AROMATICITY_RDKIT,
     ) -> Generator[tuple[tuple[Mol], tuple[Mol]], None, None]:
        '''Keep reacting and fragmenting a pair of monomers until all reactive sites have been reacted
        Returns fragment pairs at each step of the chain propagation process'''
        reactants = monomers # initialize reactive pair with monomers
        while self.rxn_schema.reactants_are_compatible(reactants):
            adducts : list[Mol] = []
            fragments : list[Mol] = []
            for adduct in self.rxn_schema.react(
                    reactants,
                    repetitions=1,
                    keep_map_labels=True, # can't clear map numbers yet, otherwise intermonomer bond finder would have nothing to work with
                    sanitize_ops=sanitize_ops,
                    aromaticity_model=aromaticity_model,
                    _suppress_reactant_validation=True, # avoid double-validation since we required it as a precheck for this protocol
                ) : 
                # DEVNOTE: consider doing fragmentation on the combined molecule made up of all products?
                for fragment in self.fragment_strategy.produce_fragments(adduct, separate=True):
                    clear_atom_props(fragment, in_place=True) # essential to avoid reaction mapping info from prior steps from contaminating future ones
                    clear_bond_props(fragment, in_place=True)
                    sanitize(fragment, sanitize_ops=sanitize_ops, aromaticity_model=aromaticity_model, in_place=True) # apply same cleanup ops to fragments as to adduct
                    fragments.append(fragment)

                if clear_map_nums: # NOTE : CRITICAL that this be done after fragmentation step, which RELIES on map numbers being present
                    clear_atom_map_nums(adduct, in_place=True)
                adducts.append(adduct)
                    
            yield tuple(adducts), tuple(fragments) # yield the adduct Mol and any subsequent resulting reactive fragments
            reactants = fragments # set fragments from current round of polymerization as reactants for next round
            
    def propagate_pooled(
            self,
            monomers : Sequence[Mol],
            n_rxn_steps_max : int=5,
            clear_map_nums : bool=True,
            sanitize_ops : SanitizeFlags=SANITIZE_ALL,
            aromaticity_model : AromaticityModel=AROMATICITY_RDKIT,
        ) -> None:
        '''
        Discovers and enumerates all possible repeat unit fragments formable from a given polymerization step reaction mechanism
        
        Propagation acts on a pool of fragments reactants (initially just the "monomers" passed in) and proceeds in rounds
        where all reactible subsets are adducted and fragmented, and any previously-unseen fragments are added to the pool
        Enumeration halts either when no new fragments have been, or the set maximum number of reaction steps is reached
        
        "Uniqueness" of a fragment is assessed by its RDKit-canonicalized SMILES representation
        '''
        ...