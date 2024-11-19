'''Tabulated references of Greek letter names, prefices, and unicode characters'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

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

_GREEK_START_IDXS = { # indices where each case of the Greek alphabet starts in Unicode
    'LOWER' : 945,
    'UPPER' : 913
}

for case, idx in _GREEK_START_IDXS.items():
    globals()[f'GREEK_{case}'] = { # add dicts to global namespace
        letter_name : chr(idx + i)
            for i, letter_name in enumerate(GREEK_LETTER_NAMES)
    }