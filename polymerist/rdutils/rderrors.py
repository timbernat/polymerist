'''Custom exceptions specific to RDKit-related functionality''' 
# TODO : consider divvying up this module down if few other modules import individual errors

class SubstructMatchFailedError(Exception):
    '''Raised when molecule graph isomorphism match does not form a cover'''
    pass

class BondOrderModificationError(Exception):
    '''Raised when an invalid RDBond modification is attempted'''
    pass