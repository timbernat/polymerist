'''For conversion of RDMols back and forth between different format encodings - often imbues a desired side effect (such as 2D-projection)'''

from abc import ABC, abstractmethod, abstractproperty
from rdkit import Chem

from .rdtypes import RDMol
from ..genutils.decorators.classmod import register_subclasses

from .mapping.bijection import bijective_atom_id_iter
from .rdprops import copy_rd_props


@register_subclasses(key_attr='TAG')
class RDConverter(ABC): # TODO : add some optional sanitization measures to ensure valid output and bijection
    '''For converting an existing RDKit Molecule to and from a particular format to gain new properties'''
    @abstractproperty
    @classmethod
    def TAG(cls):
        pass

    @abstractmethod
    def _convert(self, rdmol : RDMol) -> RDMol:
        '''Implement conversion mechanism here'''
        pass

    def convert(self, rdmol : RDMol, sanitize : bool=True) -> RDMol:
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

class SMARTSConverter(RDConverter):
    TAG = 'SMARTS'
    def _convert(self, rdmol : RDMol) -> RDMol:
        return Chem.MolFromSmarts(Chem.MolToSmarts(rdmol))

class SMILESConverter(RDConverter):
    TAG = 'SMILES'
    def _convert(self, rdmol : RDMol) -> RDMol:
        return Chem.MolFromSmiles(Chem.MolToSmiles(rdmol), sanitize=False)
    
class CXSMARTSConverter(RDConverter):
    '''Similar to SMARTSConverter but preserves the 3D structure'''
    TAG = 'CXSMARTS'
    def _convert(self, rdmol : RDMol) -> RDMol:
        return Chem.MolFromSmarts(Chem.MolToCXSmarts(rdmol))

class CXSMILESConverter(RDConverter):
    '''Similar to SMILESConverter but preserves the 3D structure'''
    TAG = 'CXSMILES'
    def _convert(self, rdmol : RDMol) -> RDMol:
        return Chem.MolFromSmiles(Chem.MolToCXSmiles(rdmol), sanitize=False)

class InChIConverter(RDConverter): # TOSELF : this does not preserve atom map num ordering (how to incorporate AuxInfo?)
    TAG = 'InChI'
    def _convert(self, rdmol : RDMol) -> RDMol:
        return Chem.AddHs(Chem.MolFromInchi(Chem.MolToInchi(rdmol), removeHs=False, sanitize=False))
    
class JSONConverter(RDConverter):
    TAG = 'JSON'
    def _convert(self, rdmol : RDMol) -> RDMol:
        return Chem.rdMolInterchange.JSONToMols(Chem.MolToJSON(rdmol))[0]
