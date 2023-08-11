'''Custom data containers with useful properties'''

from collections import defaultdict


def RecursiveDict() -> defaultdict:
    '''Returns a defaultdict which can be recursively nested indefinitely'''
    return defaultdict(lambda : RecursiveDict())