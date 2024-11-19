'''Systematic Greek numerical prefixes'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from ..combinatorics.partitions import make_change_greedy

# CONSTANT REFERENCE VALUES (adapted from https://en.wikipedia.org/wiki/List_of_polygons#Systematic_polygon_names)
GREEK_PREFIXES_UNITS = {
    1   : 'mono', # TODO : add inference rules for prefixes of numbers not tabulated here
    2   : 'di',
    3   : 'tri',
    4   : 'tetra',
    5   : 'penta',
    6   : 'hexa',
    7   : 'hepta',
    8   : 'octa',
    9   : 'ennea',
}
GREEK_PREFIXES_TEENS = { # 13-19
    10 + i : GREEK_PREFIXES_UNITS[i] + 'deca'
        for i in range(3, 9+1)
}
GREEK_PREFIXES_TENS = { # 30, 40, ..., 80, 90
    10 * i : GREEK_PREFIXES_UNITS[i] + 'conta'
        for i in range(3, 9+1)
}
GREEK_PREFIXES_SPECIAL = { # cases which do not adhere to systematic naming rules
    11 : 'hendeca',
    12 : 'dodeca',
    20 : 'icosi',
}
GREEK_PREFIXES_POW_TENS = {
    10**1 : 'deca',
    10**2 : 'hecta',
    10**3 : 'chili',
    10**4 : 'myri',
    10**6 : 'mega',
}

## combined dict
GREEK_PREFIXES = { 
    **GREEK_PREFIXES_UNITS,
    **GREEK_PREFIXES_TEENS,
    **GREEK_PREFIXES_TENS,
    **GREEK_PREFIXES_SPECIAL,
    **GREEK_PREFIXES_POW_TENS
}
GREEK_PREFIXES = { # arrange in numerical order for simpler presentation
    i : GREEK_PREFIXES[i]
        for i in sorted(GREEK_PREFIXES)
}

# DYNAMIC PREFIX GENERATION
def get_greek_prefix(n : int) -> str:
    '''Determine greek prefix for a given cardinal number'''
    is_monic = (0 < n < 10)

    tokens = []
    for i, count in make_change_greedy(n, GREEK_PREFIXES).items():
        if count != 0:
            multiplier = "" if (count == 1) else GREEK_PREFIXES.get(count)
            multiplicand = "hen" if (i == 1 and not is_monic) else GREEK_PREFIXES.get(i)

            token = multiplier + multiplicand
            tokens.append(token)
    return ''.join(tokens)