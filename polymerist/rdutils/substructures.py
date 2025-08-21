'''Utilities for counting and querying molecules substructures'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from rdkit.Chem import rdqueries, Mol


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
