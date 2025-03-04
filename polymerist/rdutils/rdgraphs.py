'''Utilities for interfacing between RDKit Mols and their graph representations'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Iterable

from rdkit import Chem
import networkx as nx

from .rdprops import detailed_rdobj_info
from ..genutils.textual.casing import camel_case_to_snake_case


# CONVERTING Mols TO Graphs
DEFAULT_ATOM_PROPS = (
    'AtomMapNum',
    'AtomicNum',
    'Degree',
    'ExplicitValence',
    'FormalCharge',
    'ImplicitValence',
    'IsAromatic',
    'Isotope',
    'Mass',
    'NoImplicit',
    'NumExplicitHs',
    'NumImplicitHs',
    'NumRadicalElectrons',
    'Symbol',
    'TotalDegree',
    'TotalNumHs',
    'TotalValence',
)

DEFAULT_BOND_PROPS = (
    'BondDir',
    'BondType',
    'BondTypeAsDouble',
    'Idx',
    'IsAromatic',
    'IsConjugated',
    'Stereo',
)

def chemical_graph(mol : Chem.Mol) -> nx.Graph:
    '''Converts an RDKit molecule to its network form, with nodes numbered by atom indices
    Returned Graph does NOT contain any chemical information or properties of any bonds or atoms, only topological (i.e. connectivty) information'''
    return nx.Graph(
        (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
            for bond in mol.GetBonds()
    )

def rdmol_to_networkx(rdmol : Chem.Mol, atom_attrs : Iterable[str]=DEFAULT_ATOM_PROPS, bond_attrs : Iterable[str]=DEFAULT_BOND_PROPS) -> nx.Graph:
    '''Convert an RDKit Mol into a NetworkX Graph with nodes numbered according to atom IDs and bonds between corresponding atoms 
    Copies over all Atom and Bond Props that have been set.
    Can also port other atom and bond properties over to the resulting graph by specifying attributes to copy over'''
    G = nx.Graph()
    
    # Initialize atom nodes, with properties
    for atom in rdmol.GetAtoms():
        atom_attr_vals = {
            camel_case_to_snake_case(attr) : attr_val
                for attr, attr_val in detailed_rdobj_info(atom).items()
                    if attr in atom_attrs
        }
        atom_attr_vals.update(atom.GetPropsAsDict())
        G.add_node(atom.GetIdx(), **atom_attr_vals)

    # Link bonded atoms with edges, with corresponding properties
    for bond in rdmol.GetBonds():
        bond_attr_vals = {
            camel_case_to_snake_case(attr) : attr_val
                for attr, attr_val in detailed_rdobj_info(bond).items()
                    if attr in bond_attrs
        }
        bond_attr_vals.update(bond.GetPropsAsDict())
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), **bond_attr_vals)

    return G