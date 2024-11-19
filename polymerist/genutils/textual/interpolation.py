'''For inserting text into other text in a rules-based manner'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import re


def insert_into_text_periodic(text : str, period : int, insertion : str='\n') -> str:
    '''Takes a string of text and another "insertion" string and inserts it throughout the text every <period> characters'''
    return insertion.join(text[i:i+period] for i in range(0, len(text), period))

def insert_into_text_periodic_re(text : str, period : int, insertion : str='\n') -> str:
    '''Takes a string of text and another "insertion" string and inserts it throughout the text every <period> characters
    Same as insert_into_text_periodic(), but implemented with regular expressions (allows for more complicated logical extensions)'''
    SPACE_RE = re.compile(f'(?s)(.{{{period}}})') # double curly braces escape the f-string syntax (to use as literals in regex quanitifer)
    return re.sub(SPACE_RE, f'\\1{insertion}', text)