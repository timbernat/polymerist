'''Utilities for building new polymer structures; currently limited to linear polymers and PDB save format'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import logging
LOGGER = logging.getLogger(__name__)

from typing import Optional

import warnings
with warnings.catch_warnings(record=True): # suppress numerous and irritating mbuild deprecation warnings
    warnings.filterwarnings('ignore',  category=DeprecationWarning)
    import mbuild as mb
    from mbuild import Compound
    from mbuild.lib.recipes.polymer import Polymer as MBPolymer

from fractions import Fraction
from pathlib import Path
from collections import Counter

from rdkit import Chem

from .exceptions import InsufficientChainLength, PartialBlockSequence, MorphologyError
from .estimation import estimate_n_atoms_linear

from ..genutils.decorators.functional import allow_string_paths
from ..genutils.textual.substrings import unique_string, repeat_string_to_length

from ..rdutils.bonding.portlib import get_linker_ids
from ..rdutils.bonding.substitution import saturate_ports, hydrogenate_rdmol_ports

from ..mdtools.openmmtools.serialization import serialize_openmm_pdb
from ..polymers.monomers.repr import MonomerGroup
from ..polymers.monomers.specification import SANITIZE_AS_KEKULE


# CONVERSION
def mbmol_from_mono_rdmol(rdmol : Chem.Mol, resname : Optional[str]=None) -> tuple[Compound, list[int]]:
    '''
    Accepts a monomer-spec-compliant SMARTS string and returns an mbuild Compound and a list of the indices of atom ports
    If "resname" is provided, will assign that name to the mBuild Compound returned
    '''
    linker_ids = [i for i in get_linker_ids(rdmol)] # record indices of ports - MUST unpack generator for mbuild compatibility
    
    # create port-free version of molecule which RDKit can embed without errors
    prot_mol = hydrogenate_rdmol_ports(rdmol, in_place=False)
    # prot_mol = saturate_ports(rdmol) # TOSELF : custom, port-based saturation methods are not yet ready for deployment - yield issues in RDKit representation under-the-hood 
    Chem.SanitizeMol(prot_mol, sanitizeOps=SANITIZE_AS_KEKULE) # ensure Mol is valid (avoids implicitValence issues)
    
    mb_compound = mb.conversion.from_rdkit(prot_mol) # native from_rdkit() method actually appears to preserve atom ordering
    if resname is not None:
        mb_compound.name = resname

    return mb_compound, linker_ids

@allow_string_paths
def mbmol_to_openmm_pdb(
        pdb_path : Path,
        mbmol : Compound, 
        num_atom_digits : int=2,
        resname_map : Optional[dict[str, str]]=None,
    ) -> None:
    '''Save an MBuild Compound into an OpenMM-compatible PDB file'''
    if resname_map is None: # avoid mutable default
        resname_map = {'RES' : 'Pol'} 

    traj = mbmol.to_trajectory() # first convert to MDTraj representation (much more infor-rich format)
    omm_top, omm_pos = traj.top.to_openmm(), traj.openmm_positions(0) # extract OpenMM representations of trajectory

    serialize_openmm_pdb(
        pdb_path,
        topology=omm_top,
        positions=omm_pos,
        uniquify_atom_ids=True,
        num_atom_id_digits=num_atom_digits,
        resname_map=resname_map
    )
    
# TODO: deduplify PDB atom anme and residue numbering code against serialize_openmm_pdb()
def mbmol_to_rdmol(
        mbmol : Compound,
        uniquify_atom_ids : bool=False,
        num_atom_id_digits : int=2,
        resname_map : Optional[dict[str, str]]=None
    ) -> Chem.Mol:
    '''Convert an mBuild Compound into an RDKit Mol, with correct atom coordinates and PDB residue info'''
    if resname_map is None:
        resname_map = {}
    
    rdmol = mbmol.to_rdkit()
    conformer = Chem.Conformer()
    conformer.Set3D(True)

    atom_id : int = 0
    element_counter = Counter()
    for resnum, mb_monomer in enumerate(mbmol.children, start=1):
        resname = resname_map.get(mb_monomer.name, mb_monomer.name[:3]) # if no remapping is found, just take first 3 chars
        # NOTE: the order of monomers and atoms within those monomers were added in the same order as iterated over here...
        #... so the atom indices **SHOULD** be in the correct order (hate that this even might be uncertain)
        for mbatom in mb_monomer.particles(): 
            conformer.SetAtomPosition(atom_id, 10*mbatom.pos.astype(float)) # conveert from nm to angstrom

            # set PDB residue info if monomer hierarchy is present
            if mbatom != mb_monomer: # for Compounds with a flat hierarchy, the children and particles of children will coincide
                symbol = mbatom.element.symbol
                atom_ser_id = element_counter[symbol]
                atom_ser_str = f'{atom_ser_id:0{num_atom_id_digits}d}' if uniquify_atom_ids else '  ' # double space keeps column justification correct when non-unique
                atom_name = f' {symbol}{atom_ser_str}' # need a leading space to get column alignment in PDB compliant with spec
                
                pdb_info = Chem.AtomPDBResidueInfo(
                    atomName=atom_name, 
                    residueName=resname,
                    residueNumber=resnum,
                    chainId='1',
                    isHeteroAtom=True,
                )
                element_counter[symbol] += 1 # only increment AFTER prior value has been assigned to the current atom
                rdmol.GetAtomWithIdx(atom_id).SetPDBResidueInfo(pdb_info)
            
            atom_id += 1 # TODO: this is an awful waay of keeping track of atom indices, see if there's a more secure way to do this
    conf_id = rdmol.AddConformer(conformer)
    
    return rdmol

# LINEAR POLYMER BUILDING
def build_linear_polymer(
        monomers : MonomerGroup,
        n_monomers : int,
        sequence : str='A',
        allow_partial_sequences : bool=False,
        add_Hs : bool=False,
        energy_minimize : bool=False,
    ) -> MBPolymer:
    '''Accepts a dict of monomer residue names and SMARTS (as one might find in a monomer JSON)
    and a degree of polymerization (i.e. chain length in number of monomers)) and returns an mbuild Polymer object'''
    # 0) DETERMINE THE ORIENTATION AND NUMBER OF TERMINAL MONOMERS, SUPPLYING THIS IF AN INVALID DEFINITION IS PROVIDED - DEV: consider moving this logic into MonomerGroup
    if monomers.has_valid_linear_term_orient: 
        term_orient = monomers.term_orient
        LOGGER.info(f'Using pre-defined terminal group orientation {term_orient}')
    else:
        term_orient = {
            resname : orient
                for (resname, rdmol), orient in zip(monomers.iter_rdmols(term_only=True), ['head', 'tail'])
        }
        LOGGER.warning(f'No valid terminal monomer orientations defined; autogenerated orientations "{term_orient}"; USER SHOULD VERIFY THIS YIELDS A CHEMICALLY-VALID POLYMER!')

    # 1) DETERMINE NUMBER OF SEQUENCE REPEATS NEEDED TO MEET TARGET NUMBER OF MONOMER UNITS (IF POSSIBLE) - DEV: consider making a separate function
    block_size = len(sequence)
    n_mono_term = len(term_orient)           # number of terminal monomers are actually present and well-defined
    n_mono_middle = n_monomers - n_mono_term # number of terminal monomers needed to reach target; in a linear chain, all monomers are either middle or terminal
    if n_mono_middle < 0:
        raise InsufficientChainLength(f'Registered number of terminal monomers exceeds requested chain length ({n_monomers}-mer chain can\'t possibly contain {n_mono_term} terminal monomers)')
    
    n_seq_whole : int     # number of full sequence repeats to reach a number of monomers less than or equal to the target
    n_symbols_remaining : int # number of any remaining symbols in sequence (i.e. monomers) needed to close the gap to the target (allowed to be 0 if target is a multiple of the sequence length)
    n_seq_whole, n_symbols_remaining = divmod(n_mono_middle, block_size) 
    print(n_seq_whole, n_symbols_remaining)
    
    if n_symbols_remaining != 0: # a whole number of sequence repeats (including possibly 0) plus some fraction of a full block sequence
        if not allow_partial_sequences:
            raise PartialBlockSequence(
                f'Partial polymer block sequence required to meet target number of monomers ("{sequence[:n_symbols_remaining]}" prefix of sequence "{sequence}"). ' \
                'If this is acceptable, set "allow_partial_sequences=True" and try calling build routine again'
            )    
        sequence_selected = repeat_string_to_length(sequence, target_length=n_mono_middle, joiner='')
        n_seq_repeats = 1 # just repeat the entire mixed-fraction length sequence (no full sequence repeats to exploit)
        LOGGER.warning(
            f'Target number of monomers is achievable WITH a partial {n_symbols_remaining}/{block_size} sequence repeat; ' \
            f'({n_seq_whole}*{block_size} [{sequence}] + {n_symbols_remaining} [{sequence[:n_symbols_remaining]}]) middle monomers + {n_mono_term} terminal monomers = {n_monomers} total monomers'
        )
    else: # for a purely-whole number of block sequence repeats
        if n_seq_whole < 1: # NOTE: if it were up to me, this would be < 0 to allow dimers, but mBuild has forced by hand
            raise InsufficientChainLength(
                f'{n_monomers}-monomer chain cannot accomodate both {n_mono_term} end groups AND at least 1 middle monomer sequence'
            )
        sequence_selected = sequence # NOTE: rename here is for clarity, and for consistency with partial sequence case
        n_seq_repeats = n_seq_whole
        LOGGER.info(
            f'Target chain length achievable with {n_seq_repeats} whole block(s) of the sequence "{sequence_selected}"; ' \
            f'({n_seq_repeats}*{block_size} [{sequence_selected}]) middle monomers + {n_mono_term} terminal monomers = {n_monomers} total monomers'
        )
    print(sequence_selected, n_seq_repeats)
    
    # 2) REGISTERING MONOMERS TO BE USED FOR CHAIN ASSEMBLY
    monomers_selected = MonomerGroup() # used to track and estimate sized of the monomers being used for building
    ## 2A) ADD MIDDLE MONOMERS TO CHAIN
    chain = MBPolymer() 
    for (resname, middle_monomer), symbol in zip(
            monomers.iter_rdmols(term_only=False),
            unique_string(sequence_selected, preserve_order=True), # only register a new monomer for each appearance of a new indicator in the sequence
        ): # zip with sequence limits number of middle monomers to length of block sequence
        LOGGER.info(f'Registering middle monomer {resname} (block identifier "{symbol}")')
        mb_monomer, linker_ids = mbmol_from_mono_rdmol(middle_monomer, resname=resname)
        chain.add_monomer(compound=mb_monomer, indices=linker_ids)
        monomers_selected.monomers[resname] = monomers.monomers[resname]

    ## 2B) ADD TERMINAL MONOMERS TO CHAIN
    term_iters = { # need to convert to iterators to allow for generator-like advancement (required for term group selection to behave as expected)
        resname : iter(rdmol_list)   # made necessary by annoying list-bound structure of current substructure spec
            for resname, rdmol_list in monomers.rdmols(term_only=True).items() 
    }
    for resname, head_or_tail in term_orient.items():
        term_monomer = next(term_iters[resname]) # will raise StopIteration if the terminal monomer in question is empty
        LOGGER.info(f'Registering terminal monomer {resname} (orientation "{head_or_tail}")')
        mb_monomer, linker_ids = mbmol_from_mono_rdmol(term_monomer, resname=resname)
        chain.add_end_groups(compound=mb_monomer, index=linker_ids.pop(), label=head_or_tail, duplicate=False) # use single linker ID and provided head-tail orientation
        monomers_selected.monomers[resname] = monomers.monomers[resname]

    # 3) ASSEMBLE AND RETURN CHAIN
    if not monomers_selected.is_linear: # verify the selected monomers actually define a linear polymer
        raise MorphologyError('Linear polymer building does not support non-linear monomer input')
    
    n_atoms_est = estimate_n_atoms_linear(monomers_selected, n_monomers) # TODO: create new MonomerGroup with ONLY the registered monomers to guarantee accuracy
    LOGGER.info(f'Assembling linear {n_monomers}-mer chain (estimated {n_atoms_est} atoms)')
    chain.build(n_seq_repeats, sequence=sequence_selected, add_hydrogens=add_Hs) # "-2" is to account for term groups (in mbuild, "n" is the number of times to replicate just the middle monomers)
    for atom in chain.particles():
        atom.charge = 0.0 # initialize all atoms as being uncharged (gets rid of pesky blocks of warnings)
    LOGGER.info(f'Successfully assembled linear {n_monomers}-mer chain (exactly {chain.n_particles} atoms)')
    
    if energy_minimize:
        LOGGER.info('Energy-minimizing chain to find more stable conformer')
        chain.energy_minimize()
        LOGGER.info('Energy minimization completed')

    return chain