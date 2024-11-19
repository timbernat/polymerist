'''Conversion tools between various programming language cases (https://en.wikipedia.org/wiki/Letter_case#Use_within_programming_languages)'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

def snake_case_to_camel_case(varname : str) -> str:
    '''Convert a name from Snake Case to Camel Case
    E.g. name_of_a_thing -> NameOfAThing'''
    return ''.join(word.capitalize() for word in varname.split('_'))

def camel_case_to_snake_case(varname : str) -> str:
    '''Convert a name from Camel Case to Snake Case
    E.g. NameOfAThing -> name_of_a_thing'''
    cap_idxs = [i for i, char in enumerate(varname) if char.isupper()]
    return '_'.join(
        varname[i_start:i_end].lower()
            for i_start, i_end in zip(cap_idxs, cap_idxs[1:]+[None])
    ) 