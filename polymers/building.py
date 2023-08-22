'''Utilities for building new polymer structures; currently limited to linear polymers and PDB save format'''

import logging
LOGGER = logging.getLogger(__name__)

from rdkit import Chem

import mbuild as mb
from mbuild import Compound
from mbuild.lib.recipes.polymer import Polymer as MBPolymer

from .exceptions import SubstructMatchFailedError, CrosslinkingError
from .estimation import estimate_chain_len_linear
from ..monomers.repr import MonomerGroup

from ..rdutils.reactions.ports import get_port_ids, hydrogenate_rdmol_ports
from ..rdutils.rdconvert import SMILESConverter
from ..rdutils.rdtypes import RDMol


# CONVERSION
def mbmol_from_mono_rdmol(rdmol : RDMol) -> tuple[Compound, list[int]]:
    '''Accepts a monomer-spec-compliant SMARTS string and returns an mbuild Compound and a list of the indices of hydrogen ports
    Alternative implementation which bypasses the unreliable substructure match (fails on aromatics even with Kekulization)'''
    rdmol = SMILESConverter().convert(rdmol)
    Chem.SanitizeMol(rdmol)

    port_ids = get_port_ids(rdmol) # record indices of ports
    prot_mol = hydrogenate_rdmol_ports(rdmol, in_place=False) # replace ports with Hs to give complete fragments
    mb_compound = mb.conversion.from_rdkit(prot_mol) # native from_rdkit() method actually appears to preserve atom ordering

    return mb_compound, port_ids

def mbmol_from_mono_rdmol_legacy(rdmol : RDMol) -> tuple[Compound, list[int]]:
    '''Accepts a monomer-spec-compliant SMARTS string and returns an mbuild Compound and a list of the indices of hydrogen ports'''
    orig_port_ids = get_port_ids(rdmol) # record indices of ports
    prot_mol = hydrogenate_rdmol_ports(rdmol, in_place=False) # replace ports with Hs to give complete fragments
    mono_smiles = Chem.MolToSmiles(prot_mol) # NOTE : CRITICAL that this be done AFTER hydrogenation (to avoid having ports in SMILES, which mbuild doesn't know how to handle)
    
    mb_compound = mb.load(mono_smiles, smiles=True)
    mb_ordered_rdmol = Chem.MolFromSmiles(mb_compound.to_smiles()) # create another molecule which has the same atom ordering as the mbuild Compound
    mb_ordered_rdmol = Chem.AddHs(mb_ordered_rdmol) # mbuild molecules don't have explicit Hs when converting to SMILES (although luckily AddHs adds them in the same order)
    Chem.Kekulize(mb_ordered_rdmol, clearAromaticFlags=True) # need to kekulize in order for aromatic bonds to be properly substructure matched (otherwise, ringed molecules are unsupported)

    mb_isomorphism = mb_ordered_rdmol.GetSubstructMatch(prot_mol) # determine mapping between original and mbuild atom indices
    if not mb_isomorphism: # ensure that the structures were in fact able to be matched before attempting backref map
        raise SubstructMatchFailedError
    mb_port_ids = [mb_isomorphism[idx] for idx in orig_port_ids]  # find the indices of the ports in the mbuild molecule

    return mb_compound, mb_port_ids


# LINEAR POLYMER BUILDING
def build_linear_polymer(monomers : MonomerGroup, DOP : int, sequence : str='A', add_Hs : bool=False, reverse_term_labels : bool=False) -> MBPolymer:
    '''Accepts a dict of monomer residue names and SMARTS (as one might find in a monomer JSON)
    and a degree of polymerization (i.e. chain length in number of monomers)) and returns an mbuild Polymer object'''
    if not monomers.is_linear:
        raise CrosslinkingError('Linear polymer building does not support non-linear monomer input')

    chain = MBPolymer() 
    term_labels = ['head', 'tail'] # mbuild requires distinct labels in order to include both term groups
    if reverse_term_labels:
        term_labels = term_labels[::-1]

    for (resname, monomer) in monomers.iter_rdmols:
        try: # attempt both methods for mbuild conversion
            mb_monomer, port_ids = mbmol_from_mono_rdmol(monomer)
        except:
            mb_monomer, port_ids = mbmol_from_mono_rdmol_legacy(monomer) # plan to deprecate, only here to maximize structure coverage
        print(mb_monomer, port_ids)
        
        if MonomerGroup.is_terminal(monomer):
            chain.add_end_groups(compound=mb_monomer, index=port_ids[0], label=term_labels.pop(), duplicate=False)
        else:
            chain.add_monomer(compound=mb_monomer, indices=port_ids)

    LOGGER.info(f'Building linear polymer chain with {DOP} monomers ({estimate_chain_len_linear(monomers, DOP)} atoms)')
    chain.build(DOP - 2, sequence=sequence, add_hydrogens=add_Hs) # "-2" is to account for term groups (in mbuild, "n" is the number of times to replicate just the middle monomers)
    for atom in chain.particles():
        atom.charge = 0.0 # initialize all atoms as being uncharged (gets risk of pesky blocks of warnings)

    return chain