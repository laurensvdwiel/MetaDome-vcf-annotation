import logging
import inspect
from dev_settings import LOGGER_NAME

def initLogging(logging_level=logging.DEBUG, print_to_console=False):
    """Initializes the logging throughout the project. If debug mode is active,
     it will also initialize a logger for the debug stream."""
    
    # ensure we are not initializing the same handler twice
    if len(logging.getLogger(LOGGER_NAME).handlers) == 0:
        # implement the basic configuration to use to write to the log files
        logging.basicConfig(filename = LOGGER_NAME+'.log',
                            format='%(asctime)s %(levelname)-2s : %(name)s %(pathname)s %(module)s.%(funcName)s %(processName)s (ID: %(process)d) %(threadName)s (ID: %(thread)d) %(message)s',
                            filemode='w', # If option 'w' is set, it will overwrite the log file with each run
                            level = logging_level) # what is written to the log file always contains debugging lines
        
        # ensure we capture warnings within packages
        logging.captureWarnings(True)
        
        # check if we are running the code via PyDev or python debugger (PDB)
        if isDebugging() or print_to_console: 
            # define a Handler which writes INFO messages or higher to the sys.stderr
            console = logging.StreamHandler()
            
            # set the logging level
            console.setLevel(logging_level)
            
            # set a format which is simpler to interpret for console use
            formatter = logging.Formatter('%(asctime)s %(levelname)-2s %(module)s.%(funcName)s: %(message)s (%(processName)s (Thread-ID: %(thread)d))')
             
            # tell the handler to use this format
            console.setFormatter(formatter)
            
            # ensure we push all custom log messages to the console
            logging.getLogger(LOGGER_NAME).addHandler(console)
            
            # ensure we push all warning to the console
            logging.getLogger('py.warnings').addHandler(console)
            
def isDebugging():
    """See: http://stackoverflow.com/questions/333995/how-to-detect-that-python-code-is-being-executed-through-the-debugger
    This checks on run time if the code is being executed in a debugging state in either pydev or python debugger (PDB)"""
    for frame in inspect.stack():
        if frame[1].endswith("pydevd.py") or frame[1].endswith("pdb.py"):
            return True
    return False