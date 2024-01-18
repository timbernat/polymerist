'''Query chemical data trough the PubChem PUG REST API (https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest)'''

import pubchempy as pcp


PUBCHEMPY_ERRORS = tuple( # need to make tuple to allow for multi-Exception catching without further processing
    getattr(pcp, objname)
        for objname in dir(pcp)
            if 'Error' in objname
)

PUBCHEMPY_ERROR_MAP = {
    exc.__name__ : exc
        for exc in PUBCHEMPY_ERRORS
}