'''Query chemical data from the NIH NCI CACTUS Chemical Identifier Resolver (https://cactus.nci.nih.gov/)'''

import requests


# CACTUS constants
CACTUS_URL = 'https://cactus.nci.nih.gov/chemical/structure'
CACTUS_PROPERTIES = {
    'Standard InChIKey'                            : 'stdinchikey',
    'Standard InChI'                               : 'stdinchi',
    'SMILES'                                       : 'smiles',
    'FICTS Identifier'                             : 'ficts',
    'FICuS Identifier'                             : 'ficus',
    'uuuuu Identifier'                             : 'uuuuu',
    'Cactvs HASHISY'                               : 'hashisy',
    'file?format=sdf'                              : 'SD File',
    'Names'                                        : 'names',
    'IUPAC Name'                                   : 'iupac_name',
    'CAS Registry Number'                          : 'cas',
    'ChemSpider ID'                                : 'chemspider_id',
    'GIF Image'                                    : 'image',
    'TwirlyMol (3D)'                               : 'twirl',
    'Molecular Weight'                             : 'mw',
    'Chemical Formula'                             : 'formula',
    'Number of Hydrogen Bond Donors'               : 'h_bond_donor_count',
    'Number of Hydrogen Bond Acceptors'            : 'h_bond_acceptor_count',
    'Number of Hydrogen Bond Acceptors and Donors' : 'h_bond_center_count',
    'Number of Rule of 5 Violations'               : 'rule_of_5_violation_count',
    'Number of Freely Rotatable Bonds'             : 'rotor_count',
    'Number of Effectively Rotatable Bonds'        : 'effective_rotor_count',
    'Number of Rings'                              : 'ring_count',
    'Number of Ring Systems'                       : 'ringsys_count',
}

class NoCACTUSDataFound(Exception):
    '''Raised when CACTUS does not have data for an otherwise valid query'''
    pass


# UTILITY FUNCTIONS
def _get_CACTUS_prop_key(prop_query : str) -> str:
    '''Determine URL extension for a given query based, or raise Exception if invalid'''
    if prop_query in CACTUS_PROPERTIES.values():
        return prop_query
    elif prop_query in CACTUS_PROPERTIES.keys():
        return CACTUS_PROPERTIES[prop_query] # check for dropdown alias
    else:
        raise ValueError(
            f'Cannot query property "{prop_query}" from NCI Resolver'
            f'\nProperty must be one of the following literals : {list(CACTUS_PROPERTIES.values())}'
            f'\nOR one of the following value aliases : {list(CACTUS_PROPERTIES.keys())}'
        )

def query_NIH_CACTUS(smiles : str, prop : str) -> str:
    '''Format URL and send request for data to NCI Resolver'''
    prop = _get_CACTUS_prop_key(prop) # sanitize and validate query

    url = f'{CACTUS_URL}/{smiles}/{prop}'
    response = requests.get(url)
    
    try:
        response.raise_for_status()
        return response.text
    except requests.HTTPError:
        raise NoCACTUSDataFound


