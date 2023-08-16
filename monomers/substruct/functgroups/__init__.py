'''SMARTS-based queries for functional groups and other chemical signatures'''

# Creating module-specific logger
import logging
LOGGER = logging.getLogger(__name__)

# locating Path to module for SMARTS data I/O
from pathlib import Path
_MODULE_PATH = Path(__file__).parent

# load/generating functional group smarts table
from . import _daylight_scrape
import pandas as pd


_FN_GRP_TABLE_NAME : str = 'fn_group_smarts'
FN_GRP_TABLE_PATH = _MODULE_PATH / f'{_FN_GRP_TABLE_NAME}.csv'

if not FN_GRP_TABLE_PATH.exists(): # if data table is missing, scrape data back off of Daylight SMARTS sight and save
    LOGGER.warning(F'No functional group SMARTS data from LUT found on system; regenerating from {_daylight_scrape.DAYLIGHT_URL}')
    FN_GROUP_TABLE = _daylight_scrape.scrape_SMARTS()
    FN_GROUP_TABLE.sort_values('category', inplace=True)  # sort by category
    FN_GROUP_TABLE.reset_index(drop=True, inplace=True)   # renumber by new order
    FN_GROUP_TABLE.to_csv(FN_GRP_TABLE_PATH, index=False) # save for future lookup
else:
    LOGGER.info('Loading functional group SMARTS data from LUT')
    FN_GROUP_TABLE = pd.read_csv(FN_GRP_TABLE_PATH) # otherwise, read data from table on import