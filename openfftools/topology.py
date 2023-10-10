'''Utilities for serializing, converting, and extracting info from OpenFF topologies'''

from pathlib import Path
from ast import literal_eval

from openff.toolkit import ToolkitRegistry
from openff.toolkit.topology import Atom, Molecule, Topology

from . import GTR
from ..genutils.iteration import asiterable
from ..genutils.fileutils.pathutils import dotless
from ..genutils.decorators.functional import allow_string_paths, optional_in_place

from ..rdutils.rdprops import RDPROP_GETTERS, RDPROP_SETTERS
from ..rdutils.rdtypes import RDAtom, RDMol


# TOPOLOGY INFO
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


# MOLECULE METADATA SERIALIZATION / DESERIALIZATION
@optional_in_place
def package_atom_metadata(offmol : Molecule) -> None:
    '''Collect atom metadata into a serializable property of the parent Molecule'''
    packaged_mdat = {
        atom.molecule_atom_index : atom.metadata
            for atom in offmol.atoms
                if atom.metadata
    }

    if packaged_mdat: # only assign if metadata is present (cleaner than assigning empty dict)
        offmol.properties['metadata'] = packaged_mdat

@optional_in_place
def unpackage_atom_metadata(offmol : Molecule) -> None:
    '''Reassign atom metadata from a "packaged" metadata property dict belonging to the parent Molecule'''
    packaged_mdat = offmol.properties.pop('metadata', None) # lookup and remove serialized data, or None if not present
    if not packaged_mdat: # handles both the null (i.e. NoneType) case and the empty dict case
        return
    
    if isinstance(packaged_mdat, str):
        packaged_mdat = literal_eval(packaged_mdat) # de-stringify if necessary

    for atom_id, metadata in packaged_mdat.items():
        offmol.atoms[atom_id].metadata.update(metadata)


# FILE I/O FUNCTIONS
@allow_string_paths
def save_molecule(path : Path, offmol : Molecule, toolkit_registry : ToolkitRegistry=GTR) -> None: # TODO : generalize to work for bytes w/ opened file
    '''Syntactic sugar for annoying suffix re-specification when saving OpenFF Molecules'''
    offmol.to_file(str(path), file_format=dotless(path), toolkit_registry=toolkit_registry)

@allow_string_paths
def topology_to_sdf(path : Path, offtop : Topology, toolkit_registry : ToolkitRegistry=GTR) -> None:
    '''Save an OpenFF Topology to an SDF file file atom metadata preserved'''
    assert(path.suffix == '.sdf')

    with path.open('w') as sdf_file:# TODO : maybe change to append mode instead?
        for mol in offtop.molecules:
            serial_mol = package_atom_metadata(mol, in_place=False)  # packageage metadata for serialization WITHOUT disturbing the original molecule
            serial_mol.to_file(sdf_file, file_format='SDF', toolkit_registry=toolkit_registry) # NOTE : not using save_molecule() here, as that function does not currently support bytes-like input (needed for multiple mols in same file)

@allow_string_paths
def topology_from_sdf(path : Path, *args, **kwargs) -> Topology:
    '''Load an OpenFF Topology from an SDF file, assigning metadata and partial charges if stored'''
    assert(path.suffix == '.sdf')

    return Topology.from_molecules(
        unpackage_atom_metadata(mol, in_place=False) # reassign metadata if serialized
            for mol in asiterable(Molecule.from_file(path, *args, **kwargs)) # asiterable() gracefully handles the singleton case
    )