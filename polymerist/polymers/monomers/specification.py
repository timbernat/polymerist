'''Implementations of the canonical monomer substructure SMARTS specification defined in https://doi.org/10.26434/chemrxiv-2023-f2zxd-v2'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import logging
LOGGER = logging.getLogger(__name__)

import re
from typing import Union

from rdkit import Chem
from rdkit.Chem import QueryAtom

from ...smileslib.primitives import is_valid_SMILES, is_valid_SMARTS, RDKIT_QUERYBONDS_BY_BONDTYPE
from ...rdutils.labeling import molwise


# CHEMICAL INFO SPECIFICATION
SANITIZE_AS_KEKULE = (Chem.SANITIZE_ALL & ~Chem.SANITIZE_SETAROMATICITY) # sanitize everything EXCEPT reassignment of aromaticity

def expanded_SMILES(
        smiles : str,
        assign_map_nums : bool=True,
        start_from : int=1,
        kekulize : bool=True,
    ) -> str:
    '''
    Expands and clarifies the chemical information contained within a passed SMILES string
    namely explicit hydrogens and bond orders, and (optionally) kekulized aromatic bonds and atom map numbers
    '''
    assert(is_valid_SMILES(smiles))
    
    rdmol = Chem.MolFromSmiles(smiles, sanitize=True)
    rdmol = Chem.AddHs(rdmol, addCoords=True)
    if assign_map_nums:
        rdmol = molwise.assign_ordered_atom_map_nums(rdmol, start_from=start_from)
    
    if kekulize:
        Chem.Kekulize(rdmol, clearAromaticFlags=True)
    Chem.SanitizeMol(rdmol)

    return Chem.MolToSmiles(rdmol, kekuleSmiles=kekulize, allBondsExplicit=True, allHsExplicit=True)


# REGEX TEMPLATES FOR COMPLIANT SMARTS
COMPLIANT_ATOM_SMARTS = re.compile( # spec-compliant query
    r'\[' \
    r'(?P<isotope>\d?)' \
    r'#(?P<atomic_num>\d+?)' \
    r'(D(?P<degree>\d{1}))?' \
    r'([+-](?P<formal_charge>\d+))?'
    r':(?P<atom_map_num>\d+?)' \
    r'\]'
)

ABERRANT_ATOM_SMARTS = re.compile( # spec-compliant query which has been mangled by RDKit
    r'\[' \
    r'#(?P<atomic_num>\d+?)' \
    r'(&(?P<isotope>\d+?)\*)?' \
    r'(&D(?P<degree>\d{1}))?' \
    r'(&[+-](?P<formal_charge>\d+))?'
    r':(?P<atom_map_num>\d+?)' \
    r'\]'
)

def chem_info_from_match(match : re.Match) -> dict[str, Union[int, str, None]]:
    '''Generate chemical information dict (with proper types) from an atom query regex match'''
    atom_info = match.groupdict()
    for key, value in atom_info.items():
        if isinstance(value, str) and value.isdigit():
            atom_info[key] = int(value) # convert parsed strings to ints where possible

    return atom_info


# SMARTS ATOM QUERY GENERATION
def compliant_atom_query_from_info(
        atomic_num : int,
        degree : int,
        atom_map_num : int,
        formal_charge : int=0,
        isotope : int=0,
        as_atom : bool=False
    ) -> Union[str, QueryAtom]:
    '''Construct a monomer-spec compliant atom SMARTS string directly from chemical information'''
    if not isotope: # handles when isotope is literal 0 or NoneType
        isotope = "" # non-specific isotope is not explicitly written in string (left empty)
    
    if atomic_num == 0: # define dummy/wild atoms to be linkers
        # atom_query = f'[{isotope}*:{atom_map_num}]'
        atom_query = f'[{isotope}#{atomic_num}:{atom_map_num}]'
    else:
        atom_query = f'[{isotope}#{atomic_num}D{degree}{formal_charge:+}:{atom_map_num}]'

    if as_atom:
        return Chem.AtomFromSmarts(atom_query)
    return atom_query

def compliant_atom_query_from_rdatom(rdatom : Chem.Atom, as_atom : bool=False) -> Union[str, QueryAtom]:
    '''Construct a monomer-spec compliant atom SMARTS string from an RDKit Atom'''
    return compliant_atom_query_from_info(
        atomic_num   = rdatom.GetAtomicNum(),
        degree       = rdatom.GetDegree(), # counts number of active bonds
        atom_map_num = rdatom.GetAtomMapNum(), # TODO : add check for nonzero map num   
        formal_charge= rdatom.GetFormalCharge(),
        isotope      = rdatom.GetIsotope(),
        as_atom=as_atom
    )

def compliant_atom_query_from_re_match(match : re.Match) -> str:
    '''Construct a monomer-spec compliant atom SMARTS string from a RegEx string match of a compliant or aberrant atom'''
    return compliant_atom_query_from_info(**chem_info_from_match(match), as_atom=False)


# CONVERSION METHODS
## DEV: add function to check whether a given SMARTS is COMPLETELY spec-compliant
def compliant_mol_SMARTS(smarts : str) -> str:
    '''Convert generic SMARTS string into a spec-compliant one'''
    # initial checks
    assert(is_valid_SMARTS(smarts))
    rdmol = Chem.MolFromSmarts(smarts)
    assert(molwise.has_fully_mapped_atoms(rdmol))
    assert(molwise.has_uniquely_mapped_atoms(rdmol))
    # TODO : add aromaticity checks
    
    # assign query info to atoms and bonds
    for atom in rdmol.GetAtoms():
        atom_query = compliant_atom_query_from_rdatom(atom, as_atom=True)
        atom.SetQuery(atom_query)

    for bond in rdmol.GetBonds():
        bond_query = RDKIT_QUERYBONDS_BY_BONDTYPE[bond.GetBondType()]
        bond.SetQuery(bond_query)

    # sanitize away RDKit artifacts
    unsanitized_smarts = Chem.MolToSmarts(rdmol) # RDKit does some mangling of queries which needs to be corrected here
    sanitized_smarts, num_repl = re.subn(
        pattern=ABERRANT_ATOM_SMARTS,
        repl=compliant_atom_query_from_re_match,
        string=unsanitized_smarts,
        count=rdmol.GetNumAtoms() # can't possibly replace more queries than there are atoms
    )
    if num_repl > 0:
        LOGGER.debug(f'Cleaned {num_repl} SMARTS atom query aberrations introduced by RDKit')
    sanitized_smarts = sanitized_smarts.replace('#0', '*') # replace explicit atom number 0 calls with star (easier to do post-processing, as #0 is easier to implement)

    return sanitized_smarts
