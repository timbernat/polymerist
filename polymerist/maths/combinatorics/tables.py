'''Tabular and printable forms of collections of combinatorial numbers'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from .numbers import binomial_coeff

def pascal(N : int) -> str:
    '''
    Generate a printable string of the first N rows of Pascal's triangle
    
    Parameters
    ----------
    N : int
        Number of rows to print
    
    Returns
    -------
    pascal_string : str
        A string consisting of each row of Pascal's triangle,
        center-justified and separated by newline characters
        
    '''
    row_width : int = None
    pascal_rows : list[str] = []
    for n in reversed(range(N)): # generate rows from longest-to-shortest
        pascal_row = ' '.join(
            f'{round(binomial_coeff(n, k))}' # NOTE: this could in principle be generalized to other 2-parameter numbers (e.g. the Stirling numbers)
                for k in range(n + 1)
        )
        
        if row_width is None:
            row_width = len(pascal_row) # set maximum row width on first row only, which is by construction the longest
        else:
            pascal_row = pascal_row.center(row_width)
        pascal_rows.append(pascal_row)
        
    return '\n'.join(pascal_rows[::-1]) # un-reverse for final return