'''For generating linear polymer structure from monomer, sequence, and chain length information'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import logging
LOGGER = logging.getLogger(__name__)

import warnings
with warnings.catch_warnings(record=True): # suppress numerous and irritating mbuild deprecation warnings
    warnings.filterwarnings('ignore',  category=DeprecationWarning)
    from mbuild.lib.recipes.polymer import Polymer as MBPolymer

from .mbconvert import mbmol_from_mono_rdmol
from .sequencing import LinearCopolymerSequencer
from ..exceptions import MorphologyError
from ..monomers.repr import MonomerGroup
from ..estimation import estimate_n_atoms_linear
from ...genutils.textual.substrings import unique_string


def build_linear_polymer(
        monomers : MonomerGroup,
        n_monomers : int,
        sequence : str='A',
        minimize_sequence : bool=True,
        allow_partial_sequences : bool=False,
        add_Hs : bool=False,
        energy_minimize : bool=False,
    ) -> MBPolymer:
    '''Accepts a dict of monomer residue names and SMARTS (as one might find in a monomer JSON)
    and a degree of polymerization (i.e. chain length in number of monomers)) and returns an mbuild Polymer object'''
    # 1) DETERMINE NUMBER OF SEQUENCE REPEATS NEEDED TO MEET TARGET NUMBER OF MONOMER UNITS (IF POSSIBLE) - DEV: consider making a separate function
    end_groups = monomers.linear_end_groups() # cache end groups so they dont need to be recalculated when registering end groups
    end_group_names = [resname for (resname, _) in end_groups.values()]
    
    sequencer = LinearCopolymerSequencer(
        sequence_kernel=sequence,
        n_repeat_units=n_monomers,
        n_repeat_units_terminal=len(end_groups)
    )
    if minimize_sequence:
        sequencer.reduce() # identify minimal subsequences
    
    sequence_compliant, n_seq_repeats = sequencer.procrustean_alignment(allow_partial_sequences=allow_partial_sequences)
    LOGGER.info(
        f'Target chain length achievable with {sequencer.describe_tally()}, ' \
        f'namely with the sequence {sequencer.describe_order(end_group_names=end_group_names)}'
    )
    sequence_unique = unique_string(sequence_compliant, preserve_order=True) # only register a new monomer for each appearance of a new, unique symbol in the sequence
    
    # 2) REGISTERING MONOMERS TO BE USED FOR CHAIN ASSEMBLY
    polymer = MBPolymer() 
    monomers_selected = MonomerGroup() # used to track and estimate sized of the monomers being used for building
    
    ## 2A) ADD MIDDLE MONOMERS TO CHAIN
    for symbol, (resname, middle_monomer) in zip(sequence_unique, monomers.iter_rdmols(term_only=False)): # zip with sequence limits number of middle monomers to length of block sequence
        LOGGER.info(f'Registering middle monomer {resname} (block identifier "{symbol}")')
        mb_monomer, linker_ids = mbmol_from_mono_rdmol(middle_monomer, resname=resname)
        polymer.add_monomer(compound=mb_monomer, indices=linker_ids)
        monomers_selected.monomers[resname] = monomers.monomers[resname]

    ## 2B) ADD TERMINAL MONOMERS TO CHAIN
    for head_or_tail, (resname, term_monomer) in end_groups.items():
        LOGGER.info(f'Registering terminal monomer {resname} (orientation "{head_or_tail}")')
        mb_monomer, linker_ids = mbmol_from_mono_rdmol(term_monomer, resname=resname)
        polymer.add_end_groups(compound=mb_monomer, index=linker_ids.pop(), label=head_or_tail, duplicate=False) # use single linker ID and provided head-tail orientation
        monomers_selected.monomers[resname] = monomers.monomers[resname]

    # 3) ASSEMBLE AND RETURN CHAIN
    if not monomers_selected.is_linear: # verify the selected monomers actually define a linear polymer
        raise MorphologyError('Linear polymer building does not support non-linear monomer input')
    
    n_atoms_est = estimate_n_atoms_linear(monomers_selected, n_monomers) # TODO: create new MonomerGroup with ONLY the registered monomers to guarantee accuracy
    LOGGER.info(f'Assembling linear {n_monomers}-mer chain (estimated {n_atoms_est} atoms)')
    
    polymer.build(n_seq_repeats, sequence=sequence_compliant, add_hydrogens=add_Hs) # "-2" is to account for term groups (in mbuild, "n" is the number of times to replicate just the middle monomers)
    for atom in polymer.particles():
        atom.charge = 0.0 # initialize all atoms as being uncharged (gets rid of pesky blocks of warnings)
    LOGGER.info(f'Successfully assembled linear {n_monomers}-mer chain (exactly {polymer.n_particles} atoms)')
    
    # 4) OPTIONALLY, PERFORM FINAL UFF ENERGY MINIMIZATION
    if energy_minimize:
        LOGGER.info('Energy-minimizing chain to find more stable conformer')
        polymer.energy_minimize()
        LOGGER.info('Energy minimization completed')

    return polymer