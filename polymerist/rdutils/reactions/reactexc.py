'''Exceptions specific to reactions'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

class BadNumberReactants(Exception):
    '''To be raised when too many or too few Mols are provided than expected'''
    pass

class ReactantTemplateMismatch(Exception):
    '''To be raised when a provided sequence of Mols does not match ChemicalReaction Reactant Templates'''
    pass

class ProductTemplateMismatch(Exception):
    '''To be raised when a provided sequence of Mols does not match ChemicalReaction Product Templates'''
    pass

class NoIntermonomerBondsFound(Exception):
    '''To be raised when search for newly-formed inter-monoer bonds fail'''
    pass