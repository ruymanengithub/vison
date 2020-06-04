"""

These functions can be used for logging information.

.. Warning:: logger is not multiprocessing safe.

:author: Sami-Matias Niemi

:version: 0.3
"""
from __future__ import print_function
import logging
import logging.handlers
#import new
from pdb import set_trace as stop
import textwrap
import string as st


def f_text_wrapper(msg):
    """ """
    textwrap.break_on_hypens = True
    width = 80
    #wmsg = []
    if type(msg) in [str, unicode]:
        wmsg = textwrap.wrap(msg, width=width)
        #wmsg = ['%s\\' % item for item in wmsg]
    elif isinstance(msg, list):
        wmsg = []
        for item in msg:
            wmsg += f_text_wrapper(item)

    return wmsg


def _info(f):
    def wrapper(self, msg, *args, **kwargs):

        msg = f_text_wrapper(msg)

        if isinstance(msg, str):
            _msg = msg
        if isinstance(msg, list):
            _msg = ''
            for item in msg:
                _msg += '%s\n' % item
        return f(self, _msg, *args, **kwargs)
    return wrapper


class _myLogger(logging.Logger):

    @_info
    def info(self, msg, *args, **kwargs):
        """
        Log 'msg % args' with severity 'INFO'.

        To pass exception information, use the keyword argument exc_info with
        a true value, e.g.

        logger.info("Houston, we have a %s", "interesting problem", exc_info=1)
        """
        if self.isEnabledFor(logging.INFO):
            self._log(logging.INFO, msg, args, **kwargs)


def setUpLogger(log_filename, loggername='logger'):
    """
    Sets up a logger.

    :param: log_filename: name of the file to save the log.
    :param: loggername: name of the logger

    :return: logger instance
    """
    # create logger
    logger = _myLogger(loggername)
    #logger = logging.Logger(loggername)
    logger.setLevel(logging.DEBUG)
    # Add the log message handler to the logger
    handler = logging.handlers.RotatingFileHandler(log_filename)
    # maxBytes=20, backupCount=5)
    # create formatter
    #formatter = logging.Formatter('%(asctime)s - %(module)s - %(funcName)s - %(levelname)s - %(message)s')
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    # add formatter to ch
    handler.setFormatter(formatter)
    # add handler to logger
    logger.addHandler(handler)

    return logger


class SimpleLogger(object):
    """
    A simple class to create a log file or print the information on screen.
    """

    def __init__(self, filename, verbose=False):
        self.file = open(filename, 'w')
        self.verbose = verbose

    def write(self, text):
        """
        Writes text either to file or screen.
        """
        print(text, file=self.file)
        if self.verbose:
            print(text)
