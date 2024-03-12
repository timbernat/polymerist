'''Utilities for interfacing between RDKit Mols and their graph representations'''

from typing import Iterable

from rdkit import Chem
import networkx as nx

from ..genutils.importutils import compile_simple_getable_attrs


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
    'PropsAsDict',
    'QueryType',
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

def rdmol_to_networkx(rdmol : Chem.Mol, atom_attrs : Iterable[str]=DEFAULT_ATOM_PROPS, bond_attrs : Iterable[str]=DEFAULT_BOND_PROPS) -> nx.Graph:
    '''Convert an RDKit Mol into a NetworkX Graph with nodes numbered according to atom IDs and bonds between corresponding atoms 
    Copies over all Atom and Bond Props that have been set.
    Can also port other atom and bond properties over to the resulting graph by specifying attributes to copy over'''
    G = nx.Graph()
    for atom in rdmol.GetAtoms():
        atom_attr_vals = {
            attr : attr_val
                for attr, attr_val in compile_simple_getable_attrs(atom, getter_str='Get').items()
                    if attr in atom_attrs
        }
        atom_attr_vals.update(atom.GetPropsAsDict())
        G.add_node(atom.GetIdx(), **atom_attr_vals)

    for bond in rdmol.GetBonds():
        bond_attr_vals = {
            attr : attr_val
                for attr, attr_val in compile_simple_getable_attrs(bond, getter_str='Get').items()
                    if attr in bond_attrs
        }
        bond_attr_vals.update(atom.GetPropsAsDict())
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), **bond_attr_vals)

    return G