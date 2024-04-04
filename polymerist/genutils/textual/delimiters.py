'''Validation and parsing of text enclosed by opening and closing delimiters (i.e. parentheses-like behavior)'''

from typing import Generator
from functools import partial
from enum import StrEnum


# PARSING AND VALIDATION METHODS FOR GENERAL OPENING/CLOSING DELIMITERS
class DelimiterBalanceState(StrEnum):
    '''For indication the error state in a delimiter character pair closure check'''
    BALANCED = 'balanced'
    PREMATURE_CLOSURE = 'premature closure'
    UNCLOSED_OPENING  = 'unclosed opening'

def check_balanced_delimiters(string : str, start_char : str, end_char : str) -> tuple[bool, DelimiterBalanceState, int]:
    '''Tests if parenthesis-like start and end delimiters within a string are balanced
    (i.e. there are just as many closures as opening, and no closure occurs before and opening)
    
    Return whether the string is balanced, and the index where an imbalance has occurred (-1 if no imbalance or at end)'''
    stack : list[int] = []
    for i, char in enumerate(string):
        if char == start_char:
            stack.append(i)
        if char == end_char:
            try:
                stack.pop()
            except IndexError:
                return False, DelimiterBalanceState.PREMATURE_CLOSURE, i
    if stack:
        return False, DelimiterBalanceState.UNCLOSED_OPENING, stack.pop() # even if no premature closures are present, imbalance can still occur if too few closures have happened by end
    return True, DelimiterBalanceState.BALANCED, -1 # -1 indicates no error here

def parse_within_delimiters(string : str, start_char : str, end_char : str) -> Generator[str, None, None]:
    '''Generates substring contained by starting and ending delimiting characters (and depth of substring) via pushdown automaton mechanism'''
    stack : list[int] = []
    for i, char in enumerate(string):
        if char == start_char:
            stack.append(i)
        if (char == end_char) and stack:
            start_idx = stack.pop()
            yield string[start_idx+1:i], len(stack)


# BOILERPLATE FUNCTIONALITY FOR COMMON PARENTHETICAL DELIMITERS
COMMON_DELIMITERS = {
    'parentheses'     : '()',
    'square_brackets' : '[]',
    'curly_brackets'  : '{}',
}

for delimiter_name, (start_char, end_char) in COMMON_DELIMITERS.items(): # register specific cases as globally-accessible functions
    globals()[f'parse_{delimiter_name}'         ] = partial(parse_within_delimiters  , start_char=start_char, end_char=end_char)
    globals()[f'check_balanced_{delimiter_name}'] = partial(check_balanced_delimiters, start_char=start_char, end_char=end_char)

def validate_common_delimiters(string : str, pointer_char : str='^') -> None:
    '''For checking whether all braces-like delimiters (i.e. parenthesis, brackets, etc.) are fully opened and not closed early
    If an invalid delimiters are found, a ValueError is raised with a helpful error message that points to the location of the fault (and the reason)'''
    for delimiter_name, (start_char, end_char) in COMMON_DELIMITERS.items():
        is_balanced, balance_state, err_idx = check_balanced_delimiters(string, start_char=start_char, end_char=end_char)
        if not is_balanced:
            indicator_str = err_idx*' ' + pointer_char
            raise ValueError(f'Imbalanced {delimiter_name} found in string (reason: {balance_state.value})\n{string}\n{indicator_str}')
validate_braces = validate_common_delimiters # alias for convenience