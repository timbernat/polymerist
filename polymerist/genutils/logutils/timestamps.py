'''Tools for formatting and recording timestamps for events'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Union
from dataclasses import dataclass

import re
from datetime import datetime


# DEFINING DEFAULT TIME FORMATS
TIMESTAMP_LOG = '%Y-%m-%d %H:%M:%S' # timestamp format to use for logging
TIMESTAMP_FMT = '%m-%d-%Y_at_%H-%M-%S_%p' # timestamp format to use for datetimes
TIMESTAMP_RE  = re.compile(r'\d{2}-\d{2}-\d{4}_at_\d{2}-\d{2}-\d{2}_\w{2}') # regex to use when searching for TIMESTAMP_FMT

# TIMESTAMPING CLASSES
@dataclass
class Timestamp:
    '''For storing information on date processing'''
    fmt_str : str = TIMESTAMP_FMT # should be formatted such that the resulting string can be safely used in a filename (i.e. no slashes)
    regex : Union[str, re.Pattern] = TIMESTAMP_RE

    def timestamp_now(self) -> str:
        '''Return a string timestamped with the current date and time (at the time of calling)'''
        return datetime.now().strftime(self.fmt_str)

    def extract_datetime(self, timestr : str) -> datetime:
        '''De-format a string containing a timestamp and extract just the timestamp as a datetime object'''
        timestamps = re.search(self.regex, timestr) # pull out JUST the datetime formatting component
        return datetime.strptime(timestamps.group(), self.fmt_str) # convert to datetime object