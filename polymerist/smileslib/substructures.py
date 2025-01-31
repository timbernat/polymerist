'''Utilities related to handling SMARTS queries'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Generator, TypeVar
L = TypeVar('L')

from rdkit import Chem
from rdkit.Chem import rdqueries, Mol


# REFERENCE TABLES FOR SPECIAL ATOM QUERIES (https://www.rdkit.org/docs/RDKit_Book.html#mol-sdf-support-and-extensions))
_special_queries : dict[str, list[str]] = { # shorthand for special, non-element queries; 
    # nonheavy queries: the H postfix throughout is for hydrogen
    'MH' : ['MH', 'metal_or_H'],
    'XH' : ['XH', 'halogen_or_H'],
    'AH' : ['AH', 'any'],
    'QH' : ['QH', 'heteratom_or_H', 'noncarbon_or_H'],
    # heavy atom-specific queries; identifiers not ending in H only match heavy atoms
    'M' : ['M', 'metal', 'metal_heavy'],
    'X' : ['X', 'halogen', 'halogen_heavy'],
    'A' : ['A', 'heavy', 'any_heavy'], # note that these are NOT SMARTS!! (e.g. "A" means "aliphatic atom" if interpreted as SMARTS)
    'Q' : ['Q', 'heteratom', 'heteratom_heavy', 'noncarbon', 'noncarbon_heavy'],
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
def matches_from_labelled_substructs(target_mol : Mol, substruct_queries : dict[L, Mol]) -> Generator[L, None, None]:
    '''
    Takes a Mol object and a dict of labelled substructures and generates the labels of all substructures found in the Mol
    
    Parameters
    ----------
    target_mol : Chem.Mol
        The RDKit molecule object to query
    substruct_queries : dict[L, Chem.Mol]
        Dict keyed with arbitrary (hashable) labels whose values are substructure query Mol objects
        
    Returns
    -------
    matching_labels : Generator[L, None, None]
        Yields labels of all substructures found in target_mol
    '''
    for match_mol_name, match_mol in substruct_queries.items():
        if target_mol.HasSubstructMatch(match_mol):
            yield match_mol_name

def match_dict_from_labelled_substructs(target_mol : Mol, substruct_queries : dict[L, Mol]) -> dict[L, bool]:
    '''
    Takes a Mol object and a dict of labelled substructures and returns a new dict keyed by substructure
    labels whose values indicate whether the corresponding substructure was found in the MOl
    
    Parameters
    ----------
    target_mol : Chem.Mol
        The RDKit molecule object to query
    substruct_queries : dict[L, Chem.Mol]
        Dict keyed with arbitrary (hashable) labels whose values are substructure query Mol objects
        
    Returns
    -------
    matching_labels_dict : dict[L, bool]
        Dict keyed by the same labels indicating whether the corresponding substructure was found in target_mol
    '''
    return { # DEVNOTE: kept implementation separate from matches_from_labelled_substructs() to avoid coupling a provide a low-memory alternative function
        match_mol_name : target_mol.HasSubstructMatch(match_mol) 
            for match_mol_name, match_mol in substruct_queries.items()
    }
