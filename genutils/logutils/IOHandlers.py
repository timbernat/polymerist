'''Tools for simplifying logging from multiple sources'''

import logging
from logging import Logger
from traceback import format_exception

from pathlib import Path
from typing import Optional, Union

from .timestamps import Timestamp, TIMESTAMP_LOG
from datetime import datetime


# DATE AND TIME FORMATTING
LOGGER_REGISTRY = lambda : logging.root.manager.loggerDict # dict of all extant Loggers through all loaded modules - NOTE : written as lambda to allow for updating 
LOG_FORMATTER = logging.Formatter('%(asctime)s.%(msecs)03d [%(levelname)-8s:%(module)16s:line %(lineno)-4d] - %(message)s', datefmt=TIMESTAMP_LOG) # TODO : make this more generic

# FILE-STREAM HANDLING CLASSES
class MultiStreamFileHandler(logging.FileHandler):
    '''Class to simplify logging file I/O given multiple logger sources providing logging input
    Automatically reports process completion & runtime if process is successful, or detailed error traceback otherwise
    
    Can spawn child processes to have multiple nested levels of logging to many partitioned output files for layered processes'''
    def __init__(self, filename : Union[str, Path], mode : str='a', encoding : Optional[str]=None, delay : bool=False, errors : Optional[str]=None, # FileHandler base args
                 loggers : Union[Logger, list[Logger]]=None, formatter : logging.Formatter=LOG_FORMATTER, proc_name : str='Process') -> None:       # args specific to this class
        super().__init__(filename, mode, encoding, delay, errors)

        self.proc_name : str = proc_name
        self.id : int = self.__hash__()  # generate unique ID number for tracking children
        
        self.setFormatter(formatter)
        self.personal_logger : Logger = logging.getLogger(str(self.id)) # create unique logger for internal error logging
        self.personal_logger.addHandler(self)

        self.parent : Optional[MultiStreamFileHandler] = None # to track whether the current process is the child of another process
        self.children : dict[int, MultiStreamFileHandler] = {} # keep track of child process; purely for debug (under normal circumstances, children unregister themselves once their process is complete)

        self._loggers = []
        if loggers is None:
            return
    
        if isinstance(loggers, Logger): # only reachable if loggers is explicitly passed
            self.register_logger(loggers) # handle the singleton logger case
        else:
            self.register_loggers(*loggers) # handle a list of loggers

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
            completion_msg = f'{self.proc_name} completed in {datetime.now() - self._start_time}\n'
            log_level = logging.INFO

        self.propogate_msg(level=log_level, msg=completion_msg) # log message at the determined appropriate level, passing up to parents if necessary
        # if self.parent is not None:                                               
        #     self.parent.children.pop(self.id) # orphan the current handler once its process is complete
        self.unregister_loggers() # prevents multiple redundant writes within the same Python session

        return True # TOSELF : return Falsy value only if errors are unhandled

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

class MultiStreamFileHandlerFlexible(MultiStreamFileHandler):
    '''
    MSFHandler which is a bit more flexible (and lazy) when it comes to log file naming and creation
    Can either pass a log file path (as normal) or a directory into which to output the logfile
    If the latter is chosen, will create a logfile based on the specified process name, optionally adding a timestamp if desired
    '''
    def __init__(self, filedir : Optional[Path]=None, filename : Optional[Union[str, Path]]=None, mode : str='a', encoding : Optional[str]=None, delay : bool=False, errors : Optional[str]=None, # FileHandler base args
                 loggers : Union[Logger, list[Logger]]=None, formatter : logging.Formatter=LOG_FORMATTER, proc_name : str='Process', write_timestamp : bool=True, timestamp : Timestamp=Timestamp()) -> None: # new args specific to this class
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