'''Utilities related to handling SMARTS queries'''

from rdkit import Chem
from rdkit.Chem import rdqueries


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

# TODO : add "magic" Prop queries (rxns, SDFile, and mapping)