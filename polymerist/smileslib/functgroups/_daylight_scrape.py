'''Backend web-scraping to (re)build SMARTS lookup table from the Daylight SMARTS official site'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from dataclasses import dataclass

import requests
from bs4 import BeautifulSoup

import pandas as pd


@dataclass(frozen=True)
class FnGroupSMARTSEntry:
    '''For encapuslating SMARTS group info from Daylight SMARTS registry'''
    category   : str
    category_desc : str

    group_type : str
    group_name : str

    SMARTS : str
    SMARTS_desc : str

DAYLIGHT_URL = 'https://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html'

def scrape_SMARTS(url : str=DAYLIGHT_URL) -> pd.DataFrame:
    '''Scrape SMARTS strings and accompanying descriptions and categories off of Daylight SMARTS official site'''
    soup = BeautifulSoup(requests.get(url).content, 'html.parser')

    entries = set()
    for desc_list in soup.find_all('dl'):
        category = desc_list.find_previous('a')['name']
        category_desc = desc_list.find_previous('h2').text

        group_type = desc_list.find_previous('h3').text
        for desc_term in desc_list.find_all('dt')[::-1]: # deal with annoying nesting by iterating over innermost terms first (i.e. in reverse)
            term = desc_term.extract() # remove tag from tree to prevent it from occurring in duplicate in higher terms
            text = list(term.stripped_strings)
            
            if len(text) == 3:
                group_name, SMARTS, SMARTS_desc = term.stripped_strings
            elif len(text) == 2: # raised when attempting to unpack with wrong number of args when no description was provided
                group_name, SMARTS, SMARTS_desc = *term.stripped_strings, ''
            else:
                pass
                # print(len(text), any('example' in w.lower() for w in text), text)

            entry = FnGroupSMARTSEntry(category, category_desc, group_type, group_name, SMARTS, SMARTS_desc)
            entries.add(entry)

    return pd.DataFrame.from_records(
            entry.__dict__
                for entry in entries
    )