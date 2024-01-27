'''Utilities related to handling SMARTS queries'''

from typing import Generator

from rdkit import Chem
from rdkit.Chem import rdqueries
from rdkit.Chem.rdchem import Mol as RDMol


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


# SUBSTRUCTURE QUERY UTILITIES
def num_substruct_queries(target_mol : RDMol, substruct_query : RDMol, *args, **kwargs) -> int:
    '''Get just the number of matching substructure queries in a target Mol'''
    return len(target_mol.GetSubstructMatches(substruct_query, *args, **kwargs)) # default "asMols=False" is fine here for length

def matching_labels_from_substruct_dict(target_mol : RDMol, substruct_queries : dict[str, RDMol]) -> Generator[str, None, None]:
    '''Takes a target RDKit Mol and a string-keyed dict of SMARTS substructure query Mols and 
    yields ONLY the keys of substructures which are found in the target'''
    for match_mol_name, match_mol in substruct_queries.items():
        if target_mol.HasSubstructMatch(match_mol):
            yield match_mol_name

def matching_dict_from_substruct_dict(target_mol : RDMol, substruct_queries : dict[str, RDMol]) -> dict[str, bool]:
    '''Takes a target RDKit Mol and a string-keyed dict of SMARTS substructure query Mols and 
    returns a dict of bools with the same keys indicating whether each match is present'''
    return {
        match_mol_name : target_mol.HasSubstructMatch(match_mol)
            for match_mol_name, match_mol in substruct_queries.items()
    }

def matching_dict_from_substruct_dict_alt(target_mol : RDMol, substruct_queries : dict[str, RDMol]) -> dict[str, bool]:
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