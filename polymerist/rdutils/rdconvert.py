'''For conversion of RDKit Mols back and forth between different format encodings - often imbues a desired side effect (such as 2D-projection)'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from abc import ABC, abstractmethod
from rdkit import Chem
from rdkit.Chem.rdchem import Mol

from ..genutils.decorators.classmod import register_subclasses, register_abstract_class_attrs

from .labeling.bijection import bijective_atom_id_iter
from .rdprops import copy_rd_props


@register_subclasses(key_attr='TAG')
@register_abstract_class_attrs('TAG')
class RDConverter(ABC): # TODO : add some optional sanitization measures to ensure valid output and bijection
    '''For converting an existing RDKit Molecule to and from a particular format to gain new properties'''
    @abstractmethod
    def _convert(self, rdmol : Mol) -> Mol:
        '''Implement conversion mechanism here'''
        pass

    def convert(self, rdmol : Mol, sanitize : bool=True) -> Mol:
        '''Tranform RDKit Mol using the selected method'''
        newmol = self._convert(rdmol)
        if sanitize:
            Chem.SanitizeMol(newmol)

        copy_rd_props(rdmol, newmol) # transfer properties at the mol- and atom-level to preserve consistency
        for atom_id_1, atom_id_2 in bijective_atom_id_iter(rdmol, newmol):
            orig_atom = rdmol.GetAtomWithIdx(atom_id_1)
            new_atom = newmol.GetAtomWithIdx(atom_id_2)

            copy_rd_props(orig_atom, new_atom)

        return newmol

class SMARTSConverter(RDConverter, TAG='SMARTS'):
    def _convert(self, rdmol : Mol) -> Mol:
        return Chem.MolFromSmarts(Chem.MolToSmarts(rdmol))

class SMILESConverter(RDConverter, TAG='SMILES'):
    def _convert(self, rdmol : Mol) -> Mol:
        return Chem.MolFromSmiles(Chem.MolToSmiles(rdmol), sanitize=False)
    
class CXSMARTSConverter(RDConverter, TAG='CXSMARTS'):
    '''Similar to SMARTSConverter but preserves the 3D structure'''
    def _convert(self, rdmol : Mol) -> Mol:
        return Chem.MolFromSmarts(Chem.MolToCXSmarts(rdmol))

class CXSMILESConverter(RDConverter, TAG='CXSMILES'):
    '''Similar to SMILESConverter but preserves the 3D structure'''
    def _convert(self, rdmol : Mol) -> Mol:
        return Chem.MolFromSmiles(Chem.MolToCXSmiles(rdmol), sanitize=False)

class InChIConverter(RDConverter, TAG='InChI'):
    # TOSELF : this does not preserve atom map num ordering (how to incorporate AuxInfo?)
    def _convert(self, rdmol : Mol) -> Mol:
        return Chem.AddHs(Chem.MolFromInchi(Chem.MolToInchi(rdmol), removeHs=False, sanitize=False))
    
class JSONConverter(RDConverter, TAG='JSON'):
    def _convert(self, rdmol : Mol) -> Mol:
        return Chem.rdMolInterchange.JSONToMols(Chem.MolToJSON(rdmol))[0]
