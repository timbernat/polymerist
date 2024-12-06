'''
Enhanced conversions to and from mbuild Compound objects which
preserve more molecular information than the utilities provided by 
'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from ...genutils.importutils.dependencies import modules_installed, MissingPrerequisitePackage

if not modules_installed('mbuild'):
    MissingPrerequisitePackage(
        importing_package_name=__spec__.name,
        use_case='Translation between chemical representations of polymers',
        install_link='https://libraries.io/conda/openbabel',
        dependency_name='openbabel',
        dependency_name_formal='the OpenBabel chemical toolbox',
    )

from typing import Optional

from pathlib import Path
from collections import Counter

from rdkit import Chem

import warnings
with warnings.catch_warnings(record=True): # suppress numerous and irritating mbuild deprecation warnings
    warnings.filterwarnings('ignore',  category=DeprecationWarning)
    from mbuild import Compound
    from mbuild.conversion import from_rdkit
    
from ..monomers.specification import SANITIZE_AS_KEKULE
from ...genutils.decorators.functional import allow_string_paths
from ...rdutils.bonding.portlib import get_linker_ids
from ...rdutils.bonding.substitution import saturate_ports, hydrogenate_rdmol_ports
from ...mdtools.openmmtools.serialization import serialize_openmm_pdb



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
    
    mb_compound = from_rdkit(prot_mol) # native from_rdkit() method actually appears to preserve atom ordering
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