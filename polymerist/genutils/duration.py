'''Utilities for representing, converting, and formatting amounts of time'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import ClassVar, Union, TypeAlias
from dataclasses import dataclass, field

from time import time
from string import Template
from datetime import timedelta

from .typetools.categorical import _union_member_factory


# TIME CONVERSION CONSTANTS
SECONDS_PER_INTERVAL = { # hard-coded version which omits dependency on OpenMM/other unit engine
    'year'        : 31_557_600.0,
    'day'         : 86_400.0,
    'hour'        : 3_600.0,
    'minute'      : 60.0,
    'second'      : 1.0,
    'millisecond' : 1E-3,
    'microsecond' : 1E-6,
}
SECONDS_PER_INTERVAL_ORDERED = { # arrange in descending order by magnitude of conversion factor
    unit_name : factor
        for unit_name, factor in sorted(SECONDS_PER_INTERVAL.items(), key=lambda x : x[1], reverse=True)
}


# TYPING AND CONVERSION
Timeable : TypeAlias = Union[int, float, timedelta]
istimeable = _union_member_factory(Timeable, 'Timeable')

def _convert_interval_to_seconds(interval : Timeable) -> float:
    '''Takes an object interpretable as a duration in seconds and returns a float or int corresponding to that interval (in seconds)'''
    if isinstance(interval, float):
        return interval
    elif isinstance(interval, int):
        return float(interval)
    elif isinstance(interval, timedelta):
        return interval.total_seconds()
    # elif isinstance(interval, Quantity): # deprecated to avoid OpenMM requirement; may reintroduce standard unit engine for polymerist has been decided
    #     if not interval.unit.is_compatible(second):
    #         raise ValueError('Quantity must have units dimensions of time to be interpreted as an interval')
    #     return interval.in_units_of(second)._value
    else:
        raise TypeError(f'Unsupported type "{interval.__class__.__name__}" for time interval breakdown')
    

# REPRESENTING DURATIONS
class TimeTemplate(Template):
    '''Like a string Template, but which uses a percent to indicate fields (much like a date formatter)'''
    delimiter : str = '%'

@dataclass
class Duration:
    '''For representing, converting, and formatting a length of time'''
    year        : int = field(default_factory=int)
    day         : int = field(default_factory=int)
    hour        : int = field(default_factory=int)
    minute      : int = field(default_factory=int)
    second      : int = field(default_factory=int)
    millisecond : int = field(default_factory=int)
    microsecond : int = field(default_factory=int)

    # formatting constants
    _FMT_ALIASES : ClassVar[dict[str, str]] = { # hard-coded aliases for units of time for string formatting
        'year'        : 'Y',
        'day'         : 'D',
        'hour'        : 'H',
        'minute'      : 'M',
        'second'      : 'S',
        'millisecond' : 's',
        'microsecond' : 'f',

    }

    _FMT_DIGITS : ClassVar[dict[str, str]] = { # hard-coded numbers of digits for formatting
        'year'        : 3,
        'day'         : 2,
        'hour'        : 2,
        'minute'      : 2,
        'second'      : 2,
        'millisecond' : 3,
        'microsecond' : 3,

    }

    # interconversion
    @classmethod
    def from_seconds(cls, interval : Timeable) -> dict[str, int]:
        '''Takes an object interpretable as a duration in seconds and returns the breakdown by years, hours, minutes, seconds, and milliseconds'''
        time_remain = _convert_interval_to_seconds(interval)

        breakdown = {}
        for unit_name, factor in SECONDS_PER_INTERVAL_ORDERED.items():
            unit_val, time_remain = divmod(time_remain, factor) 
            breakdown[unit_name] = round(unit_val) # by definition, the quotient part must be an integer, so no precision is lost here

        return cls(**breakdown)
    
    def to_seconds(self) -> float:
        '''Convert the stored duration into a number of seconds'''
        return sum(
            unit_val * SECONDS_PER_INTERVAL[unit]
                for unit, unit_val in self.__dict__.items()
        )
    
    @property
    def total_seconds(self) -> float:
        '''conversion-to-seconds alias for convenience'''
        return self.to_seconds()
    
    # strftime-like formatting
    @property
    def _fmt_dict(self) -> dict[str, str]:
        '''Convert stored times into string format more amenable to string interpolation'''
        return {
            self._FMT_ALIASES[unit] : f'{unit_val:0{self._FMT_DIGITS[unit]}d}'
                for unit, unit_val in self.__dict__.items()
        }
    
    def format(self, fmt_str : str) -> str: # TODO : expand to cover all strftime functionality in the future
        '''Format stored duration similar to datetime.strftime (https://docs.python.org/3/library/datetime.html#datetime.datetime.strftime)'''
        template = TimeTemplate(fmt_str)
        return template.substitute(**self._fmt_dict)
    fmt = format # alias for convenience

class Timer:
    '''Simple context manager for measuring how long something takes'''
    def __init__(self) -> None:
        self.start_time = time()
        self.time_taken = None

    def __enter__(self) -> 'Timer':
        return self

    def __exit__(self, type, value, traceback) -> bool:
        self.time_taken = time() - self.start_time
        return True