'''Utilities for extracting properties from and converting OpenFF Topology objects'''

from openff.toolkit.topology import Atom, Molecule, Topology

from ..rdutils.rdprops import RDPROP_GETTERS, RDPROP_SETTERS
from ..rdutils.rdtypes import RDAtom, RDMol


def get_largest_offmol(offtop : Topology) -> Molecule:
    '''Return the largest molecule in a Topology'''
    return max(offtop.molecules, key=lambda offmol : offmol.n_atoms)

# RDKIT METADATA-PRESERVING INTERCONVERSION
def copy_atom_metadata(offatom : Atom, rdatom : RDAtom, preserve_type : bool=True) -> None:
    '''Copies all attributes from the metadata dict of an OpenFF-type Atom as Props of an RDKit-type atom'''

    for key, value in offatom.metadata.items():
        if (type(value) not in RDPROP_SETTERS) or (not preserve_type): # set as string if type is unspecified or if explicitly requested to
            rdatom.SetProp(key, str(value))
        else:
            setter = getattr(rdatom, RDPROP_SETTERS[type(value)]) # use the atom's setter for the appropriate type
            setter(key, value)

def to_rdkit_with_metadata(offmol : Molecule, preserve_type : bool=True) -> RDMol:
    '''Converts an OpenFF molecule to an RDKit molecule, preserving atomic metadata'''
    rdmol = offmol.to_rdkit()
    for i, offatom in enumerate(offmol.atoms):
        copy_atom_metadata(offatom, rdmol.GetAtomWithIdx(i), preserve_type=preserve_type)

    return rdmol