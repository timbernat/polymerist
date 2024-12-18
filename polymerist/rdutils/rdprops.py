'''For assigning, transferring, and removing properties of RDKit objects'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Any, Optional, TypeVar, TypeAlias, Union
from copy import deepcopy

from rdkit import Chem
from rdkit.Chem.rdchem import Atom, Bond, Mol, RWMol

from .labeling.bijection import bijective_atom_id_iter
from ..genutils.decorators.functional import optional_in_place
from ..genutils.typetools.categorical import _union_member_factory
from ..genutils.attrs import compile_argfree_getable_attrs


# RDKit-specific generics and type aliases
T = TypeVar('T')   # generic type for AtomProp attributes
RD = TypeVar('RD') # generic type to represent an RDKit object

RDObj : TypeAlias = Union[Atom, Bond, Mol, RWMol]
isrdobj = _union_member_factory(RDObj, 'RDObj')


# REFERENCE FOR "MAGIC" PROP KEYS AND DESCRIPTIONS IN RDKit (from https://www.rdkit.org/docs/RDKit_Book.html#magic-property-values) 
RDMOL_MAGIC_PROPS = {
    'MolFileComments'        : 'Read from/written to the comment line of CTABs.',
    'MolFileInfo'            : 'Read from/written to the info line of CTABs.',
    '_MolFileChiralFlag'     : 'Read from/written to the chiral flag of CTABs.',
    '_Name'                  : 'Read from/written to the name line of CTABs.',
    '_smilesAtomOutputOrder' : 'The order in which atoms were written to SMILES',
    '_smilesBondOutputOrder' : 'The order in which bonds were written to SMILES',
}

RDATOM_MAGIC_PROPS = {
    '_CIPCode'               : 'the CIP code (R or S) of the atom',
    '_CIPRank'               : 'the integer CIP rank of the atom',
    '_ChiralityPossible'     : 'set if an atom is a possible chiral center',
    '_MolFileRLabel'         : 'integer R group label for an atom, read from/written to CTABs.',
    '_ReactionDegreeChanged' : 'set on an atom in a product template of a reaction if its degree changes in the reaction',
    '_protected'             : 'atoms with this property set will not be considered as matching reactant queries in reactions',
    'dummyLabel'             : '(on dummy atoms) read from/written to CTABs as the atom symbol',
    'molAtomMapNumber'       : 'the atom map number for an atom, read from/written to SMILES and CTABs',
    'molfileAlias'           : 'the mol file alias for an atom (follows A tags), read from/written to CTABs',
    'molFileValue'           : 'the mol file value for an atom (follows V tags), read from/written to CTABs',
    'molFileInversionFlag'   : 'used to flag whether stereochemistry at an atom changes in a reaction, read from/written to CTABs, determined automatically from SMILES',
    'molRxnComponent'        : 'which component of a reaction an atom belongs to, read from/written to CTABs',
    'molRxnRole'             : 'which role an atom plays in a reaction (1=Reactant, 2=Product, 3=Agent), read from/written to CTABs',
    'smilesSymbol'           : 'determines the symbol that will be written to a SMILES for the atom',
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

# PROPERTY INSPECTION FUNCTIONS
def detailed_rdobj_info(rdobj : RDObj) -> dict[str, Any]:
    '''Extract all get-able info about a particular RDKit atom. Does NOT include any non-default Prop values (e.g. atomMapNumber)'''
    return compile_argfree_getable_attrs(rdobj, getter_re='Get', repl_str='')

def atom_ids_with_prop(rdmol : Mol, prop_name : str) -> list[int]:
    '''Returns list of atom IDs of atom which have a particular property assigned'''
    return [
        atom.GetIdx()
            for atom in rdmol.GetAtoms()
                if atom.HasProp(prop_name)
    ]

@optional_in_place
def annotate_atom_prop(rdmol : Mol, prop : str, prop_type : T=str, annotate_precision : Optional[int]=None) -> None:
    '''Labels the desired Prop for all atoms in a Mol which have it'''
    getter_type = RDPROP_GETTERS[prop_type]
    for atom in rdmol.GetAtoms():
        prop_val = getattr(atom, getter_type)(prop) # invoke type-appropriate getter on atom, with the name of the desired property
        
        if hasattr(prop_val, '__round__') and annotate_precision is not None: # only round on roundable objects, and only when
            prop_val = round(prop_val, annotate_precision)
        atom.SetProp('atomNote', str(prop_val)) # need to convert to string, as double is susceptible to float round display errors (shows all decimal places regardless of rounding)

def aggregate_atom_prop(rdmol : Mol, prop : str, prop_type : T=str) -> dict[int, T]:
    '''Collects the values of a given Prop across all atoms in an RDKit molecule'''
    getter_type = RDPROP_GETTERS[prop_type]
    return {
        atom_idx : getattr(rdmol.GetAtomWithIdx(atom_idx), getter_type)(prop) # invoke type-appropriate getter on atom, with the name of the desired property
            for atom_idx in atom_ids_with_prop(rdmol, prop_name=prop)
    }


# PROPERTY TRANSFER FUNCTIONS
def copy_rd_props(from_rdobj : RD, to_rdobj : RD) -> None: # NOTE : no need to incorporate typing info, as RDKit objects can correctly interpret typed strings
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


# PROPERTY REMOVAL FUNCTIONS
@optional_in_place
def clear_atom_props(rdmol : Mol) -> None:
    '''Wipe properties of all atoms in a molecule'''
    for atom in rdmol.GetAtoms():
        for prop_name in atom.GetPropNames():
            atom.ClearProp(prop_name)


# PROPERTY COMPARISON FUNCTIONS
def difference_rdmol(rdmol_1 : Mol, rdmol_2 : Mol, prop : str='PartialCharge', remove_map_nums : bool=True) -> Mol:
    '''
    Takes two RDKit Mols (presumed to have the same structure and atom map numbers) and the name of a property 
    whose partial charges are the differences betwwen the two Mols' charges (atomwise)
    
    Assumes that the property in question is numeric (i.e. can be interpreted as a float)
    '''
    diff_mol = deepcopy(rdmol_1) # duplicate first molecule as template
    all_deltas = []
    for rdatom_1_idx, rdatom_2_idx in bijective_atom_id_iter(rdmol_1, rdmol_2):
        rdatom_1 = rdmol_1.GetAtomWithIdx(rdatom_1_idx)
        rdatom_2 = rdmol_2.GetAtomWithIdx(rdatom_2_idx)
        diff_atom = diff_mol.GetAtomWithIdx(rdatom_1_idx) # same index, since it is a deep copy

        delta = rdatom_1.GetDoubleProp(prop) - rdatom_2.GetDoubleProp(prop)
        diff_atom.SetDoubleProp(f'Delta{prop}', delta)
        all_deltas.append(delta)

        diff_atom.ClearProp(prop) # reset property value from original atom copy to avoid confusion
        if remove_map_nums:
            diff_atom.ClearProp('molAtomMapNumber') # Remove atom map num for greater visual clarity when drawing

    diff_mol.ClearProp(prop) # reset property value from original mol copy to avoid confusion
    diff_mol.SetProp(f'Delta{prop}s', str(all_deltas)) # label stringified version of property list (can be de-stringified via ast.literal_eval)
    diff_mol.SetDoubleProp(f'Delta{prop}Min', min(all_deltas)) # label minimal property value for ease of reference
    diff_mol.SetDoubleProp(f'Delta{prop}Max', max(all_deltas)) # label maximal property value for ease of reference

    return diff_mol