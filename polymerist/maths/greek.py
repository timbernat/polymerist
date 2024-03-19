'''Tabulated references of Greek letter names, prefices, and unicode characters'''

GREEK_LETTER_NAMES = [ # names for greek character literals
    'alpha',
    'beta',
    'gamma',
    'delta',
    'epsilon',
    'zeta',
    'eta',
    'theta',
    'iota',
    'kappa',
    'lambda',
    'mu',
    'nu',
    'xi',
    'omicron',
    'pi',
    'rho',
    'sigma_end',
    'sigma',
    'tau',
    'upsilon',
    'phi',
    'chi',
    'psi',
    'omega'
]

GREEK_PREFIXES = { # adapted from https://www.georgehart.com/virtual-polyhedra/greek-prefixes.html 
    1   : 'mono', # TODO : add inference rules for prefixes of numbers not tabulated here
    2   : 'di',
    3   : 'tri',
    4   : 'tetra',
    5   : 'penta',
    6   : 'hexa',
    7   : 'hepta',
    8   : 'octa',
    9   : 'ennea',
    10  : 'deca',
    11  : 'hendeca',
    12  : 'dodeca',
    13  : 'triskaideca',
    14  : 'tetrakaideca',
    15  : 'pentakaideca',
    16  : 'hexakaideca',
    17  : 'heptakaideca',
    18  : 'octakaideca',
    19  : 'enneakaideca',
    20  : 'icosa',
    24  : 'icositetra',
    30  : 'triconta' ,
    40  : 'tetraconta',
    50  : 'pentaconta',
    60  : 'hexaconta',
    70  : 'heptaconta',
    80  : 'octaconta',
    90  : 'enneaconta',
    100 : 'hecto' ,
}

_greek_start_idxs = { # indices where each case of the Greek alphabet starts in Unicode
    'LOWER' : 945,
    'UPPER' : 913
}

for case, idx in _greek_start_idxs.items():
    globals()[f'GREEK_{case}'] = { # add dicts to global namespace
        letter_name : chr(idx + i)
            for i, letter_name in enumerate(GREEK_LETTER_NAMES)
    }