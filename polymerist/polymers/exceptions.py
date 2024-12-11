'''Custom Exceptions specific to Polymers and related objects''' # TODO: go through these and purge errors which are no longer relevant (ported from polysaccharide v1)

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

# CHAIN LENGTH AND SHAPE ERRORS
class InsufficientChainLength(Exception):
    '''Raised when the polymer molecule being built is too short'''
    pass

class ExcessiveChainLength(Exception):
    '''Raised when the polymer molecule being built is too long'''
    pass

class EndGroupDominatedChain(Exception):
    '''Raised to indicate there are more end groups present in a chain than are monomer possibly allowed'''

class MorphologyError(Exception):
    '''Raised when a polymer does not have the morphology (i.e. crosslinking, molecular weight, etc) an application expects'''
    pass

# COPOLYMER SEQUENCING ERRORS
class EmptyBlockSequence(Exception):
    '''Raised when a trivial sequence of copolymer block (i.e. the empty string "") is provided when no expected'''
    pass

class PartialBlockSequence(Exception):
    '''Raised when an non-whole number of copolymer blocks is needed to reach a target chain length (and is not allowed)'''
    pass

# POLYMERIZATION MISINFORMATION ERRORS
class AlreadySolvated(Exception):
    '''Raised when attempting to add solvent to a molecule which already has solvent'''
    pass

class ChargeMismatch(Exception):
    '''Raised when attempting to merge two objects which disagree on their charging status'''
    pass

class MissingStructureData(Exception):
    '''Raised when a managed molecule has no associated structure file (e.g. PDB, SDF, etc.)'''
    pass

class MissingMonomerData(Exception):
    '''Raised when no monomer fragment information is found for a Polymer'''
    pass
