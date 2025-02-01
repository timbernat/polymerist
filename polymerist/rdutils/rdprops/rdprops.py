'''Reference and utilities common to all RDKit objects'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Any, TypeAlias, Union
from rdkit.Chem.rdchem import Atom, Bond, Mol, RWMol

from ...genutils.attrs import compile_argfree_getable_attrs
from ...genutils.decorators.functional import optional_in_place
from ...genutils.typetools.categorical import _union_member_factory

RDObj : TypeAlias = Union[Atom, Bond, Mol, RWMol]
isrdobj = _union_member_factory(RDObj, 'RDObj')


# REFERENCE FOR "MAGIC" MOL PROP KEYS (https://www.rdkit.org/docs/RDKit_Book.html#romol-mol-in-python)
RDMOL_MAGIC_PROPS = {
    'MolFileComments'        : 'Read from/written to the comment line of CTABs.',
    'MolFileInfo'            : 'Read from/written to the info line of CTABs.',
    '_MolFileChiralFlag'     : 'Read from/written to the chiral flag of CTABs.',
    '_Name'                  : 'Read from/written to the name line of CTABs.',
    '_smilesAtomOutputOrder' : 'The order in which atoms were written to SMILES',
    '_smilesBondOutputOrder' : 'The order in which bonds were written to SMILES',
}

# REFERENCE TABLES FOR ENFORCING C++ TYPING THAT RDKit ENFORCES
RDPROP_GETTERS = {
    str   : 'GetProp',
    bool  : 'GetBoolProp',
    int   : 'GetIntProp',
    float : 'GetDoubleProp'
}
RDPROP_SETTERS = {
    str   : 'SetProp',
    bool  : 'SetBoolProp',
    int   : 'SetIntProp',
    float : 'SetDoubleProp'
}

# PROPERTY INSPECTION AND TRANSFER FUNCTIONS
def detailed_rdobj_info(rdobj : RDObj) -> dict[str, Any]:
    '''Extract all get-able info about a particular RDKit atom. Does NOT include any non-default Prop values (e.g. atomMapNumber)'''
    return compile_argfree_getable_attrs(rdobj, getter_re='Get', repl_str='')

def copy_rd_props(from_rdobj : RDObj, to_rdobj : RDObj) -> None: # NOTE : no need to incorporate typing info, as RDKit objects can correctly interpret typed strings
    '''For copying properties between a pair of RDKit Atoms or Mols'''
    # NOTE : avoid use of GetPropsAsDict() to avoid errors from restrictive C++ typing
    assert isrdobj(from_rdobj) and isrdobj(to_rdobj) # verify that both objects passed are RDKit objects...
    assert type(from_rdobj) == type(to_rdobj)        # ...AND that both objects are the same type of RDKit object

    for prop in from_rdobj.GetPropNames():
        to_rdobj.SetProp(prop, from_rdobj.GetProp(prop))

@optional_in_place
def assign_props_from_dict(prop_dict : dict[str, Any], rdobj : RDObj, preserve_type : bool=True) -> None:
    '''Copies all attributes from a strng-keyed dictionary of properties as Props of an RDKit object'''
    for key, value in prop_dict.items():
        if (type(value) not in RDPROP_SETTERS) or (not preserve_type): # set as string if type is unspecified or if explicitly requested to
            rdobj.SetProp(key, str(value))
        else:
            setter = getattr(rdobj, RDPROP_SETTERS[type(value)]) # use the atom's setter for the appropriate type
            setter(key, value) # pass key and value to setter method
