'''Tools for simplifying logging from multiple sources'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import logging
from logging import Logger
from traceback import format_exception

from types import ModuleType
from typing import Iterable, Optional, Union

from pathlib import Path
from datetime import datetime

from .timestamps import Timestamp, TIMESTAMP_LOG
from ..importutils.pkgiter import iter_submodules


# DATE AND TIME FORMATTING
LOG_FORMATTER = logging.Formatter('%(asctime)s.%(msecs)03d [%(levelname)-8s:%(module)16s:line %(lineno)-4d] - %(message)s', datefmt=TIMESTAMP_LOG) # TODO : make this more generic

def get_logger_registry() -> dict[str, Logger]:
    '''Gets all registered Loggers by module'''
    return logging.root.manager.loggerDict 

def get_active_loggers() -> list[Logger]:
    return [
        possible_logger
            for possible_logger in get_logger_registry().values()
                if isinstance(possible_logger, Logger) # omits PlaceHolder objects
    ]

def submodule_loggers(module : ModuleType, recursive : bool=True, blacklist : Optional[Iterable[str]]=None, sparse : bool=True) -> dict[str, Optional[logging.Logger]]:
    '''
    Produce a dict of any Logger objects present in each submodule. Can optionally generate recursively and blacklist certain modules
    
    Parameters
    ----------
    module : ModuleType
        The "root" module to begin importing from
        Represented in the Node object returned by this function
    recursive : bool, default=True
        Whether or not to recursively import modules from subpackages and add them to the tree
    blacklist : list[str] (optional), default None
        List of module names to exclude from tree building
        If provided, will exclude any modules whose names occur in this list
    sparse : bool, default=True
        Whether to only include modules which have a Logger defined (i.e. exclude all NoneType entries from returned dict)

    Returns
    -------
    logger_registry : dict[str, Optional[logging.Logger]]
        A dict keyed by module name whose values are the corresponding Logger bound to that module
    '''
    logger_registry = {}
    for module in iter_submodules(module, recursive=recursive, blacklist=blacklist):
        full_module_name = module.__name__
        module_logger = get_logger_registry().get(full_module_name, None)
        if isinstance(module, logging.PlaceHolder): 
            continue # exclude dummy Placeholder loggers

        if not (sparse and (module_logger is None)):
            logger_registry[full_module_name] = module_logger

    return logger_registry

# FILE-STREAM HANDLING CLASSES
class MultiStreamFileHandler(logging.FileHandler):
    '''Class to simplify logging file I/O given multiple logger sources providing logging input
    Automatically reports process completion & runtime if process is successful, or detailed error traceback otherwise
    
    Can spawn child processes to have multiple nested levels of logging to many partitioned output files for layered processes'''
    def __init__(self, filename : Union[str, Path], mode : str='a', encoding : Optional[str]=None, delay : bool=False, errors : Optional[str]=None, # FileHandler base args
                 loggers : Optional[Union[str, Logger, Iterable[Logger]]]='ALL', formatter : logging.Formatter=LOG_FORMATTER, proc_name : str='Process') -> None:       # args specific to this class
        super().__init__(filename, mode, encoding, delay, errors)

        self.proc_name : str = proc_name
        self.id : int = self.__hash__()  # generate unique ID number for tracking children
        
        self.setFormatter(formatter)
        self.personal_logger : Logger = logging.getLogger(str(self.id)) # create unique logger for internal error logging
        self.personal_logger.addHandler(self)

        self.parent : Optional[MultiStreamFileHandler] = None # to track whether the current process is the child of another process
        self.children : dict[int, MultiStreamFileHandler] = {} # keep track of child process; purely for debug (under normal circumstances, children unregister themselves once their process is complete)

        if loggers is None:
            return
        
        if isinstance(loggers, str):
            if loggers.lower() == 'all':
                loggers = get_active_loggers() # done inside of init to allow updating at runtime (rather than compile time)
            else:
                return
    
        self._loggers = []
        if isinstance(loggers, Logger): # only reachable if loggers is explicitly passed
            self.register_logger(loggers) # handle the singleton logger case
        else:
            self.register_loggers(*loggers) # handle a list of loggers

    def register_logger(self, logger : Logger) -> None:
        '''Add an individual Logger to the File Stream'''
        logger.addHandler(self)
        self._loggers.append(logger)

    def unregister_logger(self, logger : Logger) -> None:
        '''Remove an individual Logger from the collection of linked Loggers'''
        logger.removeHandler(self)
        logger_idx = self._loggers.index(logger)
        self._loggers.pop(logger_idx)
    
    def register_loggers(self, *loggers : list[Logger]) -> None:
        '''Record a new Logger and add the File handler to it - enables support for multiple Logger streams to a single file'''
        for logger in loggers:
            self.register_logger(logger)

    def unregister_loggers(self) -> None:
        '''Purge all currently registered Loggers'''
        for logger in self._loggers: # TOSELF : can't just call individual unregister_logger() method, since pop will change list size while iterating (causes some Loggers to be "skipped")
            logger.removeHandler(self)
        self._loggers.clear()

    def subhandler(self, *args, **kwargs) -> 'MultiStreamFileHandler':
        '''Generate a subordinate "child" process which reports back up to the spawning "parent" process once complete'''
        child = self.__class__(*args, **kwargs) # generalizes inheritance much better
        self.children[child.id] = child # register the child for reference
        child.parent = self # assert self as the child's parent (allows child to tell if any parent handlers exist above it)

        return child
    
    def propogate_msg(self, level : int, msg : str) -> None:
        '''Propogate a logged message up through the parent tree'''
        self.personal_logger.log(level=level, msg=msg) # log message at the determined appropriate level
        if self.parent is not None:                           # if the current process is the child of another...
            self.parent.propogate_msg(level=level, msg=msg)   # ...pass the completion message up the chain to the parent

    def __enter__(self) -> 'MultiStreamFileHandler':
        self._start_time = datetime.now()
        return self
    
    def __exit__(self, exc_type, exc, trace) -> bool:
        if exc_type is not None: # unexpected error
            completion_msg = ''.join(format_exception(exc_type, exc, trace)) # format error message and traceback similar to console printout
            log_level = logging.FATAL
        else: # normal completion of context block
            completion_msg = f'Process "{self.proc_name}" completed in {datetime.now() - self._start_time}\n'
            log_level = logging.INFO

        self.propogate_msg(level=log_level, msg=completion_msg) # log message at the determined appropriate level, passing up to parents if necessary
        # if self.parent is not None:                                               
        #     self.parent.children.pop(self.id) # orphan the current handler once its process is complete
        self.unregister_loggers() # prevents multiple redundant writes within the same Python session

        return True # TOSELF : return Falsy value only if errors are unhandled

class MultiStreamFileHandlerFlexible(MultiStreamFileHandler):
    '''
    MSFHandler which is a bit more flexible (and lazy) when it comes to log file naming and creation
    Can either pass a log file path (as normal) or a directory into which to output the logfile
    If the latter is chosen, will create a logfile based on the specified process name, optionally adding a timestamp if desired
    '''
    def __init__(self, filedir : Optional[Path]=None, filename : Optional[Union[str, Path]]=None, mode : str='a', encoding : Optional[str]=None, delay : bool=False, errors : Optional[str]=None, # FileHandler base args
                 loggers : Optional[Union[str, Logger, list[Logger]]]='ALL', formatter : logging.Formatter=LOG_FORMATTER, proc_name : str='Process', write_timestamp : bool=True, timestamp : Timestamp=Timestamp()) -> None: # new args specific to this class
        if filename is None:
            if filedir is None:
                raise AttributeError('Must specify either a path to log file OR an output directory in which to generate a log file')
            assert(filedir.is_dir()) # double check that a proper directory has in fact been passed
            
            # implicit "else" if no error is raised 
            filestem = proc_name.replace(' ', '_') # use process name to generate filename; change spaces to underscores for file saving
            if write_timestamp:
                filestem += f'_{timestamp.timestamp_now()}' 
            filename = filedir / f'{filestem}.log'

        super().__init__(filename, mode, encoding, delay, errors, loggers, formatter, proc_name)

# ALIASES TO MAKE CLASS NAMES LESS UNWIELDY
MSFHandler        = MultiStreamFileHandler
MSFHandlerFlex    = MultiStreamFileHandlerFlexible