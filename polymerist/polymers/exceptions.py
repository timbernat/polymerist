'''Custom Exceptions specific to Polymers and related objects''' # TODO: go through these and purge errors which are no longer relevant (ported from polysaccharide v1)


class SubstructMatchFailedError(Exception):
    '''Raised when molecule graph isomorphism match does not form a cover'''
    pass

class InsufficientChainLengthError(Exception):
    '''Raised when the polymer molecule being built is too short'''
    pass

class ExcessiveChainLengthError(Exception):
    '''Raised when the polymer molecule being built is too long'''
    pass

class MorphologyError(Exception):
    '''Raised when a polymer does not have the morphology (i.e. crosslinking, molecular weight, etc) an application expects'''
    pass

class AlreadySolvatedError(Exception):
    '''Raised when attempting to add solvent to a molecule which already has solvent'''
    pass

class ChargeMismatchError(Exception):
    '''Raised when attempting to merge two objects which disagree on their charging status'''
    pass

class NoSimulationsFoundError(Exception):
    '''Raised when attempting to load a simulation for a managed molecule when none are present'''
    pass

class MissingStructureData(Exception):
    '''Raised when a managed molecule has no associated structure file (e.g. PDB, SDF, etc.)'''
    pass

class MissingForceFieldData(Exception):
    '''Raised when a forcefield is unspecified for a Simulation or Interchange'''
    pass

class MissingMonomerData(Exception):
    '''Raised when no monomer information is found for a Polymer'''
    pass

class MissingMonomerDataUncharged(MissingMonomerData):
    '''Raised when no monomer information WITHOUT library charges is found for a Polymer'''
    pass

class MissingMonomerDataCharged(MissingMonomerData):
    '''Raised when no monomer information WITH library charges is found for a Polymer'''
    pass
