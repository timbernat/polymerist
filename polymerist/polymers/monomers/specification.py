'''Implementations of the canonical monomer substructure SMARTS specification defined in https://doi.org/10.26434/chemrxiv-2023-f2zxd-v2'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import logging
LOGGER = logging.getLogger(__name__)

import re
from typing import Union

from rdkit import Chem
from rdkit.Chem import QueryAtom

from ...smileslib.cleanup import (
    Smiles, 
    Smarts,
    is_valid_SMARTS,
    InvalidSMARTS,
    # DEVNOTE: DON'T TOUCH! - expanded_SMILES import remains here purely for 
    # backwards-compatibility, as expanded_SMILES was originally defined in this module
    expanded_SMILES, 
)
from ...smileslib.primitives import RDKIT_QUERYBONDS_BY_BONDTYPE
from ...rdutils.chemlabel import (
    has_fully_mapped_atoms,
    has_uniquely_mapped_atoms,
    assign_ordered_atom_map_nums,
    mol_to_smiles_and_atom_permutation,
)


# REGEX TEMPLATES FOR COMPLIANT SMARTS
COMPLIANT_ATOM_SMARTS = re.compile( # pattern for specification-compliant query
    r'''    
    \[                              ## start atom entry: open bracket
    (?P<isotope>\d?)                ## isotope label: (possibly no) digits
    \#(?P<atomic_num>\d+?)          ## atomic number: literal octothorpe followed by at least 1 digit (non-greedy)
    (D(?P<degree>\d{1}))?           ## degree label: literal D (for number of explicit connections in SMARTS) and exactly one digit
    ([+-](?P<formal_charge>\d+))    ## formal charge: explicit sign (+ or -) followed by at least one digit
    :(?P<atom_map_num>\d+?)         ## atom map number: colon followed by at least one digit (empty not allowed, since explicitly-mapped atoms are reuiqred by the spec)
    \]                              ## end atom entry: close bracket
    ''',
    flags=re.VERBOSE,
)

ABERRANT_RDKIT_ATOM_SMARTS = re.compile( # pattern for a specification-compliant query which has been mangled by RDKit
    r'''
    \[                              ## start atom entry: open bracket
    \#(?P<atomic_num>\d+?)          ## atomic number: literal octothorpe followed by at least 1 digit (non-greedy)
    (&(?P<isotope>\d+?)\*)?         ## isotope label (totally optional): ampersand delimiter, one or more digits, and a literal asterisk (which RDKit for some reason inserts)
    (&D(?P<degree>\d{1}))?          ## degree label  (totally optional): ampersand delimiter, literal D (for number of explicit connections in SMARTS) and exactly one digit
    (&                              ## formal charge (totally optional): ampersand delimiter, followed by...
        (?P<charge_sign>[+-])           ### a charge sign (non-optional + or -)
        (?P<charge_magnitude>\d?)       ### a charge magnitude, as some number of digits (possibly none is valid, e.g. "+" is understood to mean "+1")
    )?  
    :(?P<atom_map_num>\d+?)         ## atom map number: colon followed by at least one digit
    \]                              ## end atom entry: close bracket
    ''',
    flags=re.VERBOSE,
)

def disambiguate_formal_charge(sign : str, magnitude : str) -> int:
    '''
    Convert a (possibly implicit) sign and magnitude of a SMILES-compliant atom
    formal charge entry to an explicit signed integer value (e.g. "+" -> +1, "--" -> -2)
    
    Parameters
    ----------
    sign : str ("", "+", or "-"
        The string (possibly empty) representing the sign of the formal charge
    magnitude : str ("" or digit string)
        The string (possibly empty) representing the magnitude of the formal charge (digits only)
        
    Returns
    -------
    formal_charge : int
        The explicit signed integer value of the formal charge
        
        Will raise ValueError if sign and magnitude passed 
        cannot be coerced into the appropriate types
    '''
    if not magnitude:
        if not sign:
            sign, magnitude = '+0' # special case for when the field is totally blank (atoms considered neutral by default)
        else:
            magnitude = '1'

    return int(f'{sign}{magnitude}') # DEVNOTE: worth checking if magnitude is digit here? will be rejected by int conversion anyway if not

def chem_info_from_match(match : re.Match) -> dict[str, Union[int, str, None]]:
    '''Generate chemical information dict (with proper types) from an atom query regex match'''
    atom_info = match.groupdict()
    atom_info['formal_charge'] = disambiguate_formal_charge(
        sign=atom_info.pop('charge_sign'),           # raise Exception f ancilliary data is not present, or else look up and remove it to avoid danling keys when passed as args
        magnitude=atom_info.pop('charge_magnitude'), # raise Exception f ancilliary data is not present, or else look up and remove it to avoid danling keys when passed as args
    )
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
        # NOTE: to simplify info fetching from regex, require that dummy atoms be explicitly labelled as atomic number 0 (rather than "*")
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
def compliant_mol_SMARTS(smarts : Union[Smiles, Smarts]) -> str:
    '''
    Convert a generic SMARTS string into a monomer specification-compliant one
    For details on specification, see https://doi-org.colorado.idm.oclc.org/10.1021/acs.jcim.3c01691
    
    Parameters
    ----------
    smarts : Union[Smiles, Smarts]
        The SMARTS (or, by virtue of superset, SMILES) string to convert

    Returns
    -------
    compliant_smarts : Smarts
        The structurally-correspondent monomer specification-compliant SMARTS string
    '''
    # initialize Mol object from passed SMARTS
    if not is_valid_SMARTS(smarts):
        raise InvalidSMARTS(f'SMARTS string "{smarts}" cannot be interpreted as a valid SMARTS pattern')

    # downconvert from general SMARTS to SMILES for structure expansion (cannot be done on QueryMols)
    ## DEVNOTE: I know these SMILES/SMARTS conversions look superfluous, but they are absolutely necessary
    ## The cast from SMARTS to SMILES allows this function to digest SMARTS while also allowing structural operations
    ## (e.g. adding Hs to the intermediate mol), which otherwise wouldn't make sense on a naive SMARTS-based query Mol
    struct_smiles, atom_order = mol_to_smiles_and_atom_permutation(
        Chem.MolFromSmarts(smarts, mergeHs=False),
        allBondsExplicit=True,
        allHsExplicit=False,
        # kekuleSmiles=True,
        kekuleSmiles=False,
        canonical=False,
    )
    ## in general, RDKit scrambled the order of atoms when writing to SMILES, so we un-scramble them both to
    ## respect the user's preference for atom order AND as a necessary condition for this function to be idempotent
    struct_mol = Chem.RenumberAtoms(
        Chem.MolFromSmiles(struct_smiles, sanitize=False),
        newOrder=atom_order, 
    )
    
    # make chemical info in atoms and bonds totally explicit
    struct_mol.UpdatePropertyCache() # required to set valence for explicit Hs check without otherwise performing sanitization
    struct_mol = Chem.AddHs(struct_mol, addCoords=False)
    Chem.Kekulize(struct_mol, clearAromaticFlags=True) # force kekulization (required by spec)
    Chem.SanitizeMol(struct_mol)#, sanitizeOps=Chem.SANITIZE_ALL & ~Chem.SANITIZE_ADJUSTHS & ~Chem.SANITIZE_KEKULIZE & ~Chem.SANITIZE_SETAROMATICITY) # ensure sanitization does not undo kekulization above

    # ensure all atoms are uniquely mapped (also required by spec); don;t overwrite perfectly-valid existing mapping if present
    # DEVNOTE: this MUST be done AFTER addition of Hs, so new unmapped H atoms don't get inserted after map number assignment
    if not (has_fully_mapped_atoms(struct_mol) and has_uniquely_mapped_atoms(struct_mol)):
        assign_ordered_atom_map_nums(struct_mol, start_from=1, in_place=True)
    
    # assign query info to atoms and bonds
    query_mol = Chem.MolFromSmarts(
        Chem.MolToSmiles(struct_mol, kekuleSmiles=True, allBondsExplicit=True, allHsExplicit=True),
        mergeHs=False,
    )
    for atom in query_mol.GetAtoms():
        atom_query = compliant_atom_query_from_rdatom(atom, as_atom=True)
        atom.SetQuery(atom_query)

    for bond in query_mol.GetBonds():
        bond_query = RDKIT_QUERYBONDS_BY_BONDTYPE[bond.GetBondType()]
        bond.SetQuery(bond_query)

    unsanitized_smarts = Chem.MolToSmarts(query_mol) # RDKit does some mangling of queries which needs to be corrected here
    
    # sanitize away RDKit artifacts
    sanitized_smarts, num_repl = re.subn(
        pattern=ABERRANT_RDKIT_ATOM_SMARTS,
        repl=compliant_atom_query_from_re_match,
        string=unsanitized_smarts,
        count=query_mol.GetNumAtoms() # can't possibly replace more queries than there are atoms
    )
    if num_repl > 0:
        LOGGER.debug(f'Cleaned {num_repl} SMARTS atom query aberrations introduced by RDKit')
    sanitized_smarts = sanitized_smarts.replace('#0', '*') # replace explicit atom number 0 calls with star (easier to do post-processing, as #0 is easier to implement)

    return sanitized_smarts
