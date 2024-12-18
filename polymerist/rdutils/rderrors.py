'''Custom exceptions specific to RDKit-related functionality''' 

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

# TODO : consider divvying up this module down if few other modules import individual errors
class SubstructMatchFailedError(Exception):
    '''Raised when molecule graph isomorphism match does not form a cover'''
    pass

class BondOrderModificationError(Exception):
    '''Raised when an invalid RDKit bond modification is attempted'''
    pass