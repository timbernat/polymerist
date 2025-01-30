'''Utilities related to handling SMARTS queries'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Generator, TypeVar
T = TypeVar('T')

from rdkit import Chem
from rdkit.Chem import rdqueries, Mol


# REFERENCE TABLES FOR SPECIAL ATOM QUERIES (https://www.rdkit.org/docs/RDKit_Book.html#mol-sdf-support-and-extensions))
_special_queries : dict[str, list[str]] = { # shorthand for special, non-element queries; 
    # nonheavy queries: the H postfix throughout is for hydrogen
    'MH' : ['metal_or_H'],
    'XH' : ['halogen_or_H'],
    'AH' : ['any'],
    'QH' : ['heteratom_or_H', 'noncarbon_or_H'],
    # heavy atom-specific queries; identifiers not ending in H only match heavy atoms
    'M' : ['metal', 'metal_heavy'],
    'X' : ['halogen', 'halogen_heavy'],
    'A' : ['heavy', 'any_heavy'], # note that these are NOT SMARTS!! (e.g. "A" means "aliphatic atom" if interpreted as SMARTS)
    'Q' : ['heteratom', 'heteratom_heavy', 'noncarbon', 'noncarbon_heavy'],
}

SPECIAL_QUERY_SMARTS : dict[str, str] = {}
SPECIAL_QUERY_ATOMS  : dict[str, Chem.QueryAtom] = {}
SPECIAL_QUERY_MOLS   : dict[str, Chem.Mol] = {}
for query_identifier, aliases in _special_queries.items():
    query_atom = getattr(rdqueries, f'{query_identifier}AtomQueryAtom')()
    query_smarts = query_atom.GetSmarts()
    query_mol = Chem.MolFromSmarts(query_smarts)
    
    for alias in aliases:
        SPECIAL_QUERY_SMARTS[alias] = query_smarts
        SPECIAL_QUERY_ATOMS[ alias] = query_atom
        SPECIAL_QUERY_MOLS[  alias] = query_mol


# COUNTING SUBSTRUCTURE QUERIES
def num_substruct_queries(target_mol : Mol, substruct_query : Mol, *args, **kwargs) -> int:
    '''Get the number of RDKit substruct matches to a SMARTS query within a given target Mol'''
    return len(target_mol.GetSubstructMatches(substruct_query, *args, **kwargs)) # default "asMols=False" is fine here for length

def num_automorphisms(substruct_query : Mol, *args, **kwargs) -> int:
    '''Get the matches a substructure query has to itself; provides measure of the degree of symmetry of the query'''
    kwargs['uniquify'] = False # force non-unique match to see how many possible ways the atoms in the substruct can be mapped to themselves
    return num_substruct_queries(substruct_query, substruct_query, *args, **kwargs)

def num_substruct_queries_distinct(target_mol : Mol, substruct_query : Mol) -> int: # TODO : also include stereochemical symmetries (via "useChirality" flag) 
    '''Get the number of distinct, non-overlapping RDKit substruct matches to a SMARTS query within a given target Mol
    Accounts for automorphic symmetries of the query to avoid double-counting distinct groups'''
    num_distinct : float = num_substruct_queries(target_mol, substruct_query, uniquify=False) # !CRITICAL! that matches be non-unique
    num_distinct /= num_automorphisms(substruct_query) # quotient out number of substructure automorphisms
    if num_distinct.is_integer(): # should be guaranteed to be an integer by Lagrange's Theorem, but good to double check here
        return int(num_distinct)
    else:
        raise ValueError('Automorphism normalization returned a non-integer number of query matches')


# MAPPING SUBSTRUCTURE QUERIES
def matching_labels_from_substruct_dict(target_mol : Mol, substruct_queries : dict[T, Mol]) -> Generator[T, None, None]:
    '''Takes a target RDKit Mol and a string-keyed dict of SMARTS substructure query Mols and 
    yields ONLY the keys of substructures which are found in the target'''
    for match_mol_name, match_mol in substruct_queries.items():
        if target_mol.HasSubstructMatch(match_mol):
            yield match_mol_name

def matching_dict_from_substruct_dict(target_mol : Mol, substruct_queries : dict[T, Mol]) -> dict[T, bool]:
    '''Takes a target RDKit Mol and a string-keyed dict of SMARTS substructure query Mols and 
    returns a dict of bools with the same keys indicating whether each match is present'''
    return {
        match_mol_name : target_mol.HasSubstructMatch(match_mol)
            for match_mol_name, match_mol in substruct_queries.items()
    }

def matching_dict_from_substruct_dict_alt(target_mol : Mol, substruct_queries : dict[T, Mol]) -> dict[T, bool]:
    '''Takes a target RDKit Mol and a string-keyed dict of SMARTS substructure query Mols and 
    returns a dict of bools with the same keys indicating whether each match is present
    
    Alternate implementation with lazy evaluation'''
    matching_dict = {
        match_mol_name : False # set False as sentinel value
            for match_mol_name in substruct_queries.keys()
    }

    for label in matching_labels_from_substruct_dict(target_mol, substruct_queries=substruct_queries):
        matching_dict[label] = True

    return matching_dict