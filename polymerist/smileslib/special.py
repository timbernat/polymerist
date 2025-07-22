'''Reference for "special" SMARTS queries not found on the periodic table'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from rdkit.Chem import rdqueries, MolFromSmarts
from rdkit.Chem.rdchem import Mol, QueryAtom


# Covers queries from https://www.rdkit.org/docs/RDKit_Book.html#mol-sdf-support-and-extensions,
# which can be implemented based on the table here: https://www.rdkit.org/docs/RDKit_Book.html#smarts-reference
_special_queries : dict[str, list[str]] = { # shorthand for special, non-element queries; 
    # nonheavy queries: the H postfix throughout is for hydrogen
    'MH' : ['MH', 'metal_or_H'],
    'XH' : ['XH', 'halogen_or_H'],
    'AH' : ['AH', 'any'],
    'QH' : ['QH', 'heteratom_or_H', 'noncarbon_or_H'],
    # heavy atom-specific queries; identifiers not ending in H only match heavy atoms
    'M'  : ['M', 'metal', 'metal_heavy'],
    'X'  : ['X', 'halogen', 'halogen_heavy'],
    'A'  : ['A', 'heavy', 'any_heavy'], # note that these are NOT SMARTS!! (e.g. "A" means "aliphatic atom" if interpreted as SMARTS)
    'Q'  : ['Q', 'heteratom', 'heteratom_heavy', 'noncarbon', 'noncarbon_heavy'],
}

SPECIAL_QUERY_SMARTS : dict[str, str] = {}
SPECIAL_QUERY_ATOMS  : dict[str, QueryAtom] = {}
SPECIAL_QUERY_MOLS   : dict[str, Mol] = {}
for query_identifier, aliases in _special_queries.items():
    query_atom = getattr(rdqueries, f'{query_identifier}AtomQueryAtom')()
    query_smarts = query_atom.GetSmarts()
    query_mol = MolFromSmarts(query_smarts)
    
    for alias in aliases:
        SPECIAL_QUERY_SMARTS[alias] = query_smarts
        SPECIAL_QUERY_ATOMS[ alias] = query_atom
        SPECIAL_QUERY_MOLS[  alias] = query_mol
