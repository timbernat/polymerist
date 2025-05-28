'''Classes for implementing reactions with respect to some set of reactant RDKit Mols'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import logging
LOGGER = logging.getLogger(__name__)

from typing import Generator, Iterable, Sequence
from dataclasses import dataclass, field

from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdmolops import SanitizeMol, SanitizeFlags, SANITIZE_ALL
from rdkit.Chem.rdmolops import AromaticityModel, AROMATICITY_RDKIT, AROMATICITY_MDL

from .reactions import AnnotatedReaction
from .fragment import IBIS, ReseparateRGroups, CutMinimumCostBondsStrategy

from ..rdprops.atomprops import clear_atom_props
from ..rdprops.bondprops import clear_bond_props
from ..chemlabel import clear_atom_map_nums, clear_atom_isotopes
from ..sanitization import sanitize_mol

from ...smileslib.cleanup import canonical_SMILES_from_mol, Smiles


@dataclass
class PolymerizationReactor:
    '''Reactor which exhaustively generates monomers fragments according to a given a polymerization mechanism'''
    rxn_schema : AnnotatedReaction
    fragment_strategy : IBIS = field(default_factory=CutMinimumCostBondsStrategy)
    
    def propagate_pooled(
            self,
            monomers : Iterable[Mol],
            rxn_depth_max : int=5,
            allow_resampling : bool=False,
            clear_map_labels : bool=True,
            clear_dummy_labels : bool=False,
            sanitize_ops : SanitizeFlags=SANITIZE_ALL,
            aromaticity_model : AromaticityModel=AROMATICITY_MDL,
        ) -> dict[Smiles, Mol]:
        '''
        Discovers and enumerates all possible repeat unit fragments formable from a given polymerization step reaction mechanism
        
        Propagation acts on a pool of fragments reactants (initially just the "monomers" passed in) and proceeds in rounds
        where all reactible subsets are adducted and fragmented, and any previously-unseen fragments are added to the pool
        Enumeration halts either when no new fragments have been found or the set maximum number of reaction step(s) is reached
        
        "Uniqueness" of a fragment is assessed by its RDKit-canonicalized SMILES representation
        '''
        reactant_pool : dict[Smiles, Mol] = { # initialize reactants as just the monomers provided
            canonical_SMILES_from_mol(monomer) : monomer
                for monomer in monomers
        } # TODO: add option for removing initial reactants from pool at end of process
        
        rxn_depth : int = 1
        while (rxn_depth <= rxn_depth_max):
            LOGGER.info(f'Enumerating fragments formable in {rxn_depth} reaction step(s) or fewer:')
            n_new_frags_found : int = 0 # this is reassigned to reset the count to 0 each round
            for reactants in self.rxn_schema.enumerate_valid_reactant_orderings(
                reactant_pool=reactant_pool, # DEVNOTE: some room for optimization here down the line in avoiding checking seen combinations up-front
                labeling_method=canonical_SMILES_from_mol, # NOTE: reactant indices would be enough for uniqueness, but chemical labels allow for tracing identity between steps
                as_mols=True, # explicitly WANT to return canonical SMILES labels, not molecule objects  
                allow_resampling=allow_resampling,
                deterministic=True,
            ):
                if reactants is None:
                    break # NOTE: this ONLY exists the reactants loop, NOT the encapsulating while loop
                
                for adduct in self.rxn_schema.react(
                    reactants,
                    repetitions=1,
                    keep_map_labels=True, # can't clear map numbers yet, otherwise intermonomer bond finder would have nothing to work with
                    sanitize_ops=sanitize_ops,
                    aromaticity_model=aromaticity_model,
                    _suppress_reactant_validation=True, # avoid double-validation since we required it as a precheck for this protocol
                ):
                    for fragment in self.fragment_strategy.produce_fragments(adduct, separate=True):
                        # clear atom and bondproperties - essential to avoid reaction mapping info from prior steps from contaminating future ones
                        clear_atom_props(fragment, in_place=True) 
                        clear_bond_props(fragment, in_place=True)
                        
                        # apply same cleanup ops to fragments as to adduct
                        sanitize_mol(fragment, sanitize_ops=sanitize_ops, aromaticity_model=aromaticity_model, in_place=True) 
                        if clear_map_labels: # NOTE : CRITICAL that this be done after fragmentation step, which RELIES on map numbers being present
                            clear_atom_map_nums(fragment, in_place=True)
                                
                        # NOTE : !ESSENTIAL! that isotope removal not be done in-place, as this to preserve dummy labels on fragments 
                        # while still ensuring chemically-identical molecules with distinct linker labels are still considered one-and-the-same
                        fragment_dummy_free = clear_atom_isotopes(fragment, in_place=False) # canonical SMILES requires 
                        canon_smi = canonical_SMILES_from_mol(fragment_dummy_free)
                        if canon_smi not in reactant_pool:
                            LOGGER.debug(f'Discovered new fragment with canonical SMILES "{canon_smi}"')
                            reactant_pool[canon_smi] = fragment_dummy_free if clear_dummy_labels else fragment # choose which labelled fragment Mol to return based on dummy label preference
                            n_new_frags_found += 1
            
            # summarizing current rxn round and checking halting conditions
            LOGGER.info(f'Found {n_new_frags_found} new fragments formable after at least {rxn_depth} reaction step(s)')
            if n_new_frags_found == 0:
                LOGGER.info(f'HALTING NORMALLY: No new reaction fragments discovered requiring {rxn_depth} reaction step(s) or more')
                break
            rxn_depth += 1
            
        else: # only called if fragment loop is not broken out of (i.e. when while condition on rxn depth becomes Falsy)
            LOGGER.warning(f'HALTING PREMATURELY: reached the configured reaction search depth limit of {rxn_depth_max} reaction step(s)')
        
        return reactant_pool
    
    def propagate(
        self,
        monomers : Sequence[Mol],
        clear_map_labels : bool=True,
        sanitize_ops : SanitizeFlags=SANITIZE_ALL,
        aromaticity_model : AromaticityModel=AROMATICITY_MDL,
     ) -> Generator[tuple[tuple[Mol, ...], tuple[Mol, ...]], None, None]:
        '''Keep reacting and fragmenting a pair of monomers until all reactive sites have been reacted
        Returns fragment pairs at each step of the chain propagation process'''
        LOGGER.warning('PolymerizationReactor.propagate() is slated for deprecation, please migrate to PolymerizationReactor.propagate_pooled() instead')
        LOGGER.warning('Note that propagate() remains totally functional, but it does not in all cases enumerate EVERY possible repeat unit fragment.')
        LOGGER.warning('Therefore, it is recommended to use propagate_pooled() for more comprehensive fragment enumeration.')
        
        reactants = monomers # initialize reactive pair with monomers
        while (reactants := self.rxn_schema.valid_reactant_ordering(reactants)) is not None: # check for and enforce compatible reactant ordering
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
                    sanitize_mol(fragment, sanitize_ops=sanitize_ops, aromaticity_model=aromaticity_model, in_place=True) # apply same cleanup ops to fragments as to adduct
                    fragments.append(fragment)

                if clear_map_labels: # NOTE : CRITICAL that this be done after fragmentation step, which RELIES on map numbers being present
                    clear_atom_map_nums(adduct)
                adducts.append(adduct)
                    
            yield tuple(adducts), tuple(fragments) # yield the adduct Mol and any subsequent resulting reactive fragments
            reactants = fragments # set fragments from current round of polymerization as reactants for next round
            