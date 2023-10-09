'''Utilities for building new polymer structures; currently limited to linear polymers and PDB save format'''

import logging
LOGGER = logging.getLogger(__name__)

from collections import Counter
from pathlib import Path

from rdkit import Chem
from openmm.app import PDBFile

import mbuild as mb
from mbuild import Compound
from mbuild.lib.recipes.polymer import Polymer as MBPolymer

from .exceptions import MorphologyError
from .estimation import estimate_chain_len_linear
from ..monomers.repr import MonomerGroup

from ..genutils.decorators.functional import allow_string_paths
from ..rdutils.amalgamation.portlib import get_linker_ids
from ..rdutils.amalgamation.bonding import saturate_ports, hydrogenate_rdmol_ports
from ..rdutils.rdtypes import RDMol


# CONVERSION
def mbmol_from_mono_rdmol(rdmol : RDMol) -> tuple[Compound, list[int]]:
    '''Accepts a monomer-spec-compliant SMARTS string and returns an mbuild Compound and a list of the indices of atom ports'''
    linker_ids = [i for i in get_linker_ids(rdmol)] # record indices of ports - MUST unpack generator for mbuild compatibility
    
    # create chemically-complete 
    prot_mol = hydrogenate_rdmol_ports(rdmol, in_place=False)
    # prot_mol = saturate_ports(rdmol) # TOSELF : custom, port-based saturation methods are not yet ready for deployment - yield issues in RDKit representation under-the-hood 
    Chem.SanitizeMol(prot_mol) # ensure Mol is valid (avoids implicitValence issues)
    mb_compound = mb.conversion.from_rdkit(prot_mol) # native from_rdkit() method actually appears to preserve atom ordering

    return mb_compound, linker_ids

@allow_string_paths
def mbmol_to_openmm_pdb(pdb_path : Path, mbmol : Compound, num_atom_digits : int=2, res_repl : dict[str, str]=None) -> None:
    '''Save an MBuild Compound into an OpenMM-compatible PDB file'''
    if res_repl is None: # avoid mutable default
        res_repl = {'RES' : 'Pol'} 

    traj = mbmol.to_trajectory() # first convert to MDTraj representation (much more infor-rich format)
    omm_top, omm_pos = traj.top.to_openmm(), traj.openmm_positions(0) # extract OpenMM representations of trajectory

    # modifiy atom info to make readable
    counter = Counter() # for keeping track of the running index of each distinct element
    for atom in omm_top.atoms():
        symbol = atom.element.symbol
        idx = counter[symbol]
        atom.name = f'{symbol}{idx:0{num_atom_digits}d}' # extend atom name with ordered integer with specified number of digits (including leading zeros)
        counter[symbol] += 1

        repl_res_name = res_repl.get(atom.residue.name, None) # lookup current residue name to see if a replacement is called for
        if repl_res_name is not None:
            atom.residue.name = repl_res_name

    with pdb_path.open('w') as file:
        PDBFile.writeFile(omm_top, omm_pos, file)


# LINEAR POLYMER BUILDING
def build_linear_polymer(monomers : MonomerGroup, DOP : int, term_orient : dict[str, str], sequence : str='A', add_Hs : bool=False) -> MBPolymer:
    '''Accepts a dict of monomer residue names and SMARTS (as one might find in a monomer JSON)
    and a degree of polymerization (i.e. chain length in number of monomers)) and returns an mbuild Polymer object'''
    if not monomers.is_linear:
        raise MorphologyError('Linear polymer building does not support non-linear monomer input')
    
    if (term_values := sorted(term_orient.values())) != ['head', 'tail']: # sorting removes order-dependence
        raise ValueError(f'Invalid dict values {term_values}; must provide term group orientation as dict with residue names as keys and "head" and "tail" as values')

    chain = MBPolymer() 
    for (resname, monomer) in monomers.iter_rdmols:
        mb_monomer, linker_ids = mbmol_from_mono_rdmol(monomer)
        if MonomerGroup.is_terminal(monomer):
            chain.add_end_groups(compound=mb_monomer, index=linker_ids.pop(), label=term_orient[resname], duplicate=False) # use single linker ID and provided head-tail orientation
        else:
            chain.add_monomer(compound=mb_monomer, indices=linker_ids)

    LOGGER.info(f'Building linear polymer chain with {DOP} monomers ({estimate_chain_len_linear(monomers, DOP)} atoms)')
    chain.build(DOP - 2, sequence=sequence, add_hydrogens=add_Hs) # "-2" is to account for term groups (in mbuild, "n" is the number of times to replicate just the middle monomers)
    for atom in chain.particles():
        atom.charge = 0.0 # initialize all atoms as being uncharged (gets risk of pesky blocks of warnings)

    return chain