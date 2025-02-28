'''For automating RDKit molecule sanitization and chemical property cleanup'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Callable, Concatenate, Generator, Iterable, Optional, ParamSpec, TypeAlias, TypeVar, Union
P = ParamSpec('P')
from functools import wraps

from rdkit import Chem
from rdkit.Chem import Mol
from rdkit.Chem.rdmolops import SanitizeMol, SanitizeFlags, SANITIZE_ALL, SANITIZE_NONE
from rdkit.Chem.rdmolops import AromaticityModel


def sanitize_mol_outputs(mol_func : Callable[P, Union[Mol, Iterable[Mol]]]) -> Callable[Concatenate[Union[None, int, SanitizeFlags], Optional[AromaticityModel], P], Union[Mol, tuple[Mol, ...]]]:
    '''
    Decorator which injects molecule sanitization capability into a function which returns RDKit Mols
    Acts on functions which are assumed to return a single mol object, or an iterable containing Mols
    Will
    
    Decorated function will return a single mol in the former case, or a tuple of mols in the latter
    
    By default, performs NO sanitization; all sanitization operations must be explicitly specified!
    '''
    @wraps(mol_func)
    def wrapped_func(
            *args : P.args,
            sanitize_ops : Union[None, int, SanitizeFlags]=SANITIZE_NONE,
            aromaticity_model : Optional[AromaticityModel]=None,
            **kwargs : P.kwargs,
        ) -> Union[Mol, tuple[Mol, ...]]:
        # ensure sanitize operation specification is valid for target inputs
        if sanitize_ops is None: # enforce correct typing from deliberate (or accidental) NoneType pass
            sanitize_ops = SANITIZE_NONE
        
        if aromaticity_model is not None: # unset since sanitization will always be the LAST step before returning (don't want to undo aromaticity model setting)
            sanitize_ops = sanitize_ops & ~SanitizeFlags.SANITIZE_KEKULIZE & ~SanitizeFlags.SANITIZE_SETAROMATICITY
            
        # determine if return is singular, iterable, or invalid
        outputs = mol_func(*args, **kwargs)
        if isinstance(outputs, Mol):
            is_singular = True
            outputs = (outputs,)
        elif isinstance(outputs, Iterable):
            is_singular = False
            outputs = tuple(outputs) # cast to tuple in advance; all changes from here are made in-place
        else:
            raise TypeError(f'Expected wrapped function to return etiher a Chem.Mol object or an iterable of Chem.Mol, not {type(outputs)}')
        
        # perform sanitization
        for mol in outputs:
            if not isinstance(mol, Mol):
                continue # for now, tolerates other types of objects in output stream
            
            if aromaticity_model is not None:
                Chem.Kekulize(mol, clearAromaticFlags=True)
                Chem.SetAromaticity(mol, model=aromaticity_model)
            Chem.SanitizeMol(mol, sanitizeOps=sanitize_ops)
        return outputs[0] if is_singular else outputs
    return wrapped_func