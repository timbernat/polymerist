'''Utilities for building new polymer structures; currently limited to linear polymers and PDB save format'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import logging
LOGGER = logging.getLogger(__name__)

import warnings
with warnings.catch_warnings(record=True): # suppress numerous and irritating mbuild deprecation warnings
    warnings.filterwarnings('ignore',  category=DeprecationWarning)
    import mbuild as mb
    from mbuild import Compound
    from mbuild.lib.recipes.polymer import Polymer as MBPolymer

from pathlib import Path
from rdkit import Chem

from .exceptions import InsufficientChainLengthError, MorphologyError
from .estimation import estimate_n_atoms_linear

from ..genutils.decorators.functional import allow_string_paths
from ..genutils.textual.substrings import uniquify_str

from ..rdutils.bonding.portlib import get_linker_ids
from ..rdutils.bonding.substitution import saturate_ports, hydrogenate_rdmol_ports

from ..mdtools.openmmtools.serialization import serialize_openmm_pdb
from ..polymers.monomers.repr import MonomerGroup
from ..polymers.monomers.specification import SANITIZE_AS_KEKULE


# CONVERSION
def mbmol_from_mono_rdmol(rdmol : Chem.Mol) -> tuple[Compound, list[int]]:
    '''Accepts a monomer-spec-compliant SMARTS string and returns an mbuild Compound and a list of the indices of atom ports'''
    linker_ids = [i for i in get_linker_ids(rdmol)] # record indices of ports - MUST unpack generator for mbuild compatibility
    
    # create port-free version of molecule which RDKit can embed without errors
    prot_mol = hydrogenate_rdmol_ports(rdmol, in_place=False)
    # prot_mol = saturate_ports(rdmol) # TOSELF : custom, port-based saturation methods are not yet ready for deployment - yield issues in RDKit representation under-the-hood 
    Chem.SanitizeMol(prot_mol, sanitizeOps=SANITIZE_AS_KEKULE) # ensure Mol is valid (avoids implicitValence issues)
    mb_compound = mb.conversion.from_rdkit(prot_mol) # native from_rdkit() method actually appears to preserve atom ordering

    return mb_compound, linker_ids

@allow_string_paths
def mbmol_to_openmm_pdb(pdb_path : Path, mbmol : Compound, num_atom_digits : int=2, res_repl : dict[str, str]=None) -> None:
    '''Save an MBuild Compound into an OpenMM-compatible PDB file'''
    if res_repl is None: # avoid mutable default
        res_repl = {'RES' : 'Pol'} 

    traj = mbmol.to_trajectory() # first convert to MDTraj representation (much more infor-rich format)
    omm_top, omm_pos = traj.top.to_openmm(), traj.openmm_positions(0) # extract OpenMM representations of trajectory

    serialize_openmm_pdb(
        pdb_path,
        topology=omm_top,
        positions=omm_pos,
        uniquify_atom_ids=True,
        num_atom_id_digits=num_atom_digits,
        resname_repl=res_repl
    )


# LINEAR POLYMER BUILDING
def build_linear_polymer(monomers : MonomerGroup, n_monomers : int, sequence : str='A', add_Hs : bool=False, energy_minimize : bool=False) -> MBPolymer:
    '''Accepts a dict of monomer residue names and SMARTS (as one might find in a monomer JSON)
    and a degree of polymerization (i.e. chain length in number of monomers)) and returns an mbuild Polymer object'''
    # 0) VERIFY THAT CHAIN ACTUAL CAN DEFINE LINEAR POLYMER
    if not monomers.is_linear:
        raise MorphologyError('Linear polymer building does not support non-linear monomer input')
    
    if monomers.has_valid_linear_term_orient: # DEV: consider moving this logic into MonomerGroup
        term_orient = monomers.term_orient
        LOGGER.info(f'Using pre-defined terminal group orientation {term_orient}')
    else:
        term_orient = {
            resname : orient
                for (resname, rdmol), orient in zip(monomers.iter_rdmols(term_only=True), ['head', 'tail']) # will raise StopIteration if fewer
        }
        LOGGER.warning(f'No valid terminal monomer orientations defined; autogenerated orientations "{term_orient}"; USER SHOULD VERIFY THIS YIELDS A CHEMICALLY-VALID POLYMER!')

    # 1) DETERMINE NUMBER OF SEQUENCE REPEATS NEEDED TO MEET TARGET NUMBER OF MONOMER UNITS (IF POSSIBLE)
    n_terminal = len(term_orient) # determine how many terminal monomers are actually present and well-defined
    block_size = len(sequence)
    
    if ((n_monomers - n_terminal) % block_size) != 0:
        raise ValueError(f'Cannot build a(n) {n_monomers}-monomer chain from any number of {block_size}-monomer blocks and {n_terminal} end groups')
    # NOTE: not explicitly forcing n_seq_reps to catch lingering float input / inexact division errors
    n_seq_reps = (n_monomers - n_terminal) // block_size # number of times to repeat the block sequence between end groups to reach the target chain length
    if n_seq_reps < 1: # NOTE: if it were up to me, this would be < 0 to allow dimers, but mBuild has forced by hand
        raise InsufficientChainLengthError(f'{n_monomers}-monomer chain has few total monomers to accomodate {n_terminal} end groups AND at least 1 middle monomer sequence')
    # TODO: consider adding support for fractional sequence lengths IFF that fraction is a rational number whose denominator divides the sequence length...
    # ...for example, could allow 5/2 * 'BACA' to be interpreted as 'BACA|BACA|BA'; 5/3 * 'BACA' would still be invalid though
    LOGGER.info(f'Target chain length achievable with {n_seq_reps} block sequence repeat(s) ({n_seq_reps}*{block_size} [{sequence}] middle monomers + {n_terminal} terminal monomers = {n_monomers} total monomers)')

    # 2) REGISTERING MONOMERS TO BE USED FOR CHAIN ASSEMBLY
    monomers_used = MonomerGroup() # used to track and estimate sized of the monomers being used
    
    ## 2A) ADD MIDDLE MONOMERS TO CHAIN
    chain = MBPolymer() 
    for (resname, middle_monomer), sequence_key in zip(
            monomers.iter_rdmols(term_only=False),
            uniquify_str(sequence, preserve_order=True), # only register a new monomer for each appearance of a new indicator in the sequence
        ): # zip with sequence limits number of middle monomers to length of block sequence
        LOGGER.info(f'Registering middle monomer {resname} (block identifier "{sequence_key}")')
        mb_monomer, linker_ids = mbmol_from_mono_rdmol(middle_monomer)
        chain.add_monomer(compound=mb_monomer, indices=linker_ids)
        monomers_used.monomers[resname] = monomers.monomers[resname]

    ## 2B) ADD TERMINAL MONOMERS TO CHAIN
    term_iters = { # need to convert to iterators to allow for generator-like advancement (required for term group selection to behave as expected)
        resname : iter(rdmol_list)   # made necessary by annoying list-bound structure of current substructure spec
            for resname, rdmol_list in monomers.rdmols(term_only=True).items() 
    }
    for resname, head_or_tail in term_orient.items():
        LOGGER.info(f'Registering terminal monomer {resname} (orientation "{head_or_tail}")')
        term_monomer = next(term_iters[resname])
        mb_monomer, linker_ids = mbmol_from_mono_rdmol(term_monomer)
        chain.add_end_groups(compound=mb_monomer, index=linker_ids.pop(), label=head_or_tail, duplicate=False) # use single linker ID and provided head-tail orientation
        monomers_used.monomers[resname] = monomers.monomers[resname]

    # 3) ASSEMBLE AND RETURN CHAIN
    n_atoms_est = estimate_n_atoms_linear(monomers_used, n_monomers) # TODO: create new MonomerGroup with ONLY the registered monomers to guarantee accuracy
    LOGGER.info(f'Assembling linear {n_monomers}-mer chain (estimated {n_atoms_est} atoms)')
    chain.build(n_seq_reps, sequence=sequence, add_hydrogens=add_Hs) # "-2" is to account for term groups (in mbuild, "n" is the number of times to replicate just the middle monomers)
    for atom in chain.particles():
        atom.charge = 0.0 # initialize all atoms as being uncharged (gets rid of pesky blocks of warnings)
    LOGGER.info(f'Successfully assembled linear {n_monomers}-mer chain (exactly {chain.n_particles} atoms)')
    
    if energy_minimize:
        LOGGER.info('Energy-minimizing chain to find more stable conformer')
        chain.energy_minimize()
        LOGGER.info('Energy minimization completed')

    return chain