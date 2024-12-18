'''Utilities related to handling SMARTS queries'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Generator, TypeVar
T = TypeVar('T')

from rdkit import Chem
from rdkit.Chem import rdqueries, Mol


# REFERENCE TABLES FOR SPECIAL ATOM TYPES
_special = { # shorthand for special, non-element queries
    'metal'           : 'M',
    'halogen'         : 'X',
    'heavy'           : 'A',
    'heavy_noncarbon' : 'Q',
}

SPECIAL_SMARTS = { 
    query_label : getattr(rdqueries, f'{symbol}AtomQueryAtom')().GetSmarts()
        for query_label, symbol in _special.items()
}

SPECIAL_QUERY_MOLS = { # TODO : make these lambda-like so that a unique object is returned on access
    query_name : Chem.MolFromSmarts(smarts)
        for query_name, smarts in SPECIAL_SMARTS.items()
}


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