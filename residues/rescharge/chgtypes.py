'''Custom typehints and classes for residue charging'''

from dataclasses import dataclass, field
from rdkit.Chem import Mol as RDMol
from ...genutils.fileutils.jsonio import JSONifiable, JSONSerializable

# from openff.toolkit import ForceField
# from openff.toolkit.typing.engines.smirnoff.parameters import LibraryChargeHandler


ChargeMap = dict[int, float] # maps substructure IDs to partial charge values

@dataclass
class ChargedResidue:
    '''Dataclass for more conveniently storing averaged charges for a residue group'''
    charges : ChargeMap
    residue_name : str
    SMARTS : str
    mol_fragment : RDMol

@dataclass
class ChargesByResidue(JSONifiable): # make serializable to JSON
    '''Class for storing substructure charge maps by residue'''
    charges : dict[str, ChargeMap] = field(default_factory=dict)

# def write_lib_chgs_from_mono_data(monomer_group : MonomerInfo, offxml_src : Path, output_path : Path) -> tuple[ForceField, list[LibraryChargeHandler]]: # TODO - refactor to accept MonomerInfo instance
#     '''Takes a monomer JSON file (must contain charges!) and a force field XML file and appends Library Charges based on the specified monomers. Outputs to specified output_path'''
#     LOGGER.warning('Generating new forcefield XML with added Library Charges')
#     assert(output_path.suffix == '.offxml') # ensure output path is pointing to correct file type
#     assert(monomer_group.has_charges) # ensure charge entries are present

#     forcefield = ForceField(offxml_src) # simpler to add library charges through forcefield API than to directly write to xml
#     lc_handler = forcefield["LibraryCharges"]

#     lib_chgs = [] #  all library charges generated from the averaged charges for each residue
#     for resname, charge_dict in monomer_group.charges.items(): # ensures no uncharged structure are written as library charges (may be a subset of the monomers structures in the file)
#         # NOTE : original implementation deprecated due to imcompatibility with numbered ports, kept in comments here for backward compatibility and debug reasons
#         # lc_entry = { # stringify charges into form usable for library charges
#         #     f'charge{cid}' : f'{charge} * elementary_charge' 
#         #         for cid, charge in charge_dict.items()
#         # } 
#         # lc_entry['smirks'] = monomer_group['monomers'][resname] # add SMIRKS string to library charge entry to allow for correct labelling
        
#         lc_entry = {}
#         rdmol = Chem.MolFromSmarts(monomer_group.monomers[resname])

#         new_atom_id = 1 # counter for remapping atom ids - NOTE : cannot start at 0, since that would denote an invalid atom
#         for atom in sorted(rdmol.GetAtoms(), key=lambda atom : atom.GetAtomMapNum()): # renumber according to map number order, NOT arbitrary RDKit atom ordering
#             if atom.GetAtomicNum(): # if the atom is not wild type or invalid
#                 old_map_num = atom.GetAtomMapNum() # TOSELF : order of operations in this clause is highly important (leave as is if refactoring!)
#                 lc_entry[f'charge{new_atom_id}'] = f'{charge_dict[old_map_num]} * elementary_charge'

#                 atom.SetAtomMapNum(new_atom_id)
#                 new_atom_id += 1; # increment valid atom index
#             else:
#                 atom.SetAtomMapNum(0) # blank out invalid atoms in SMARTS numbering

#         lc_entry['smirks'] = Chem.MolToSmarts(rdmol)  # convert renumbered mol back to SMARTS to use for SMIRNOFF charge labelling
#         lc_params = LibraryChargeHandler.LibraryChargeType(allow_cosmetic_attributes=True, **lc_entry) # must enable cosmetic params for general kwarg passing
        
#         lc_handler.add_parameter(parameter=lc_params)
#         lib_chgs.append(lc_params)  # record library charges for reference
    
#     forcefield.to_file(output_path) # write modified library charges to new xml (avoid overwrites in case of mistakes)
#     return forcefield, lib_chgs

