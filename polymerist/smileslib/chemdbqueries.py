'''For querying chemical databases for information about molecules specified by SMILES string and other structures'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import logging
LOGGER  = logging.getLogger(__name__)

from typing import Optional

from ..genutils.importutils.dependencies import modules_installed, MissingPrerequisitePackage
if not modules_installed('cirpy'):
    raise MissingPrerequisitePackage(
        importing_package_name=__spec__.name,
        use_case='Querying the NIH CACTUS Chemical Identifier Resolver (CIR)',
        install_link='https://cirpy.readthedocs.io/en/latest/guide/install.html',
        dependency_name='cirpy',
        dependency_name_formal='CIRpy',
    )
import cirpy

from ..genutils.importutils.dependencies import modules_installed, MissingPrerequisitePackage
if not modules_installed('pubchempy'):
    raise MissingPrerequisitePackage(
        importing_package_name=__spec__.name,
        use_case='Querying the PubChem Compound database',
        install_link='https://pubchempy.readthedocs.io/en/latest/guide/install.html',
        dependency_name='pubchempy',
        dependency_name_formal='PubChemPy',
    )
import pubchempy as pcp


def get_property_from_smiles(smiles : str, prop_name : str='iupac_name') -> Optional[str]: # TODO: abstract each query method via Strategy pattern
    '''Takes the SMILES string representing a molecule and attempts to fetch its IUPAC name from NIH CACTUS and/or PubChem
    Returns the fetched IUPAC name as a str, or NoneType if both queries fail'''
    # Open with NIH query (fastest method), return name if found...
    LOGGER.debug(f'Attempting query of property "{prop_name}" from NIH CACTUS')
    iupac_name = cirpy.resolve(smiles, prop_name)
    if iupac_name is not None:
        if isinstance(iupac_name, list):
            return iupac_name.pop(0)
        return iupac_name 
    
    # ...otherwise, search through PubChem Compound queries for a matching results
    pc_prop_name = pcp.PROPERTY_MAP.get(prop_name, prop_name)
    LOGGER.debug(f'Attempting query of property "{pc_prop_name}" from PubChem PUGREST')
    for prop_query in pcp.get_properties(pc_prop_name, smiles, namespace='smiles'):
        if pc_prop_name in prop_query:
            return prop_query[pc_prop_name]
    else:
        return None 