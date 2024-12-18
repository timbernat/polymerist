'''SMARTS-based queries for functional groups and other chemical signatures'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

# Creating module-specific logger
import logging
LOGGER = logging.getLogger(__name__)

# locating Path to module for SMARTS data I/O
from pathlib import Path
_MODULE_PATH = Path(__file__).parent

import pandas as pd

from ._daylight_scrape import scrape_SMARTS, FnGroupSMARTSEntry


_fgtab_name : str = 'fn_group_smarts'
_fgtab_path = _MODULE_PATH / f'{_fgtab_name}.csv'

if not _fgtab_path.exists(): # if data table is missing, scrape data back off of Daylight SMARTS sight and save
    LOGGER.warning(F'No functional group SMARTS data from LUT found on system; regenerating from {_daylight_scrape.DAYLIGHT_URL}')
    FN_GROUP_TABLE = _daylight_scrape.scrape_SMARTS()
    FN_GROUP_TABLE.sort_values('category', inplace=True)  # sort by category
    FN_GROUP_TABLE.reset_index(drop=True, inplace=True)   # renumber by new order
    FN_GROUP_TABLE.to_csv(_fgtab_path, index=False) # save for future lookup
else:
    LOGGER.info('Loading functional group SMARTS data from LUT')
    FN_GROUP_TABLE = pd.read_csv(_fgtab_path) # otherwise, read data from table on import

FN_GROUP_ENTRIES = [ # collate into list of easily queryable entries
    FnGroupSMARTSEntry(**fgdict)
        for fgdict in FN_GROUP_TABLE.to_dict(orient='records')
]