'''For conversion of RDMols back and forth between different format encodings - often imbues a desired side effect (such as 2D-projection)'''

from abc import ABC, abstractmethod, abstractproperty
from rdkit import Chem

from .rdtypes import RDMol
from ..genutils.decorators.classmod import register_subclasses


@register_subclasses(key_attr='TAG')
class RDConverter(ABC): # TODO : add some optional sanitization measures to ensure valid output and bijection
    '''For converting an existing RDKit Molecule to and from a particular format to gain new properties'''
    @abstractproperty
    @classmethod
    def TAG(cls):
        pass

    @abstractmethod
    def convert(self, rdmol : RDMol) -> RDMol:
        pass

class SMARTSConverter(RDConverter):
    TAG = 'SMARTS'
    def convert(self, rdmol : RDMol) -> RDMol:
        return Chem.MolFromSmarts(Chem.MolToSmarts(rdmol))

class SMILESConverter(RDConverter):
    TAG = 'SMILES'
    def convert(self, rdmol : RDMol) -> RDMol:
        return Chem.MolFromSmiles(Chem.MolToSmiles(rdmol), sanitize=False)
    
class CXSMARTSConverter(RDConverter):
    '''Similar to SMARTSConverter but preserves the 3D structure'''
    TAG = 'CXSMARTS'
    def convert(self, rdmol : RDMol) -> RDMol:
        return Chem.MolFromSmarts(Chem.MolToCXSmarts(rdmol))

class CXSMILESConverter(RDConverter):
    '''Similar to SMILESConverter but preserves the 3D structure'''
    TAG = 'CXSMILES'
    def convert(self, rdmol : RDMol) -> RDMol:
        return Chem.MolFromSmiles(Chem.MolToCXSmiles(rdmol), sanitize=False)

class InChIConverter(RDConverter): # TOSELF : this does not preserve atom map num ordering (how to incorporate AuxInfo?)
    TAG = 'InChI'
    def convert(self, rdmol : RDMol) -> RDMol:
        return Chem.AddHs(Chem.MolFromInchi(Chem.MolToInchi(rdmol), removeHs=False, sanitize=False))
    
class JSONConverter(RDConverter):
    TAG = 'JSON'
    def convert(self, rdmol : RDMol) -> RDMol:
        return Chem.rdMolInterchange.JSONToMols(Chem.MolToJSON(rdmol))[0]
