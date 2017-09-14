


import matplotlib
matplotlib.use("TkAgg")


from pipe import FlatFielding
from pipe.master import Pipe
from support.report import Report
from eyegore.eyegore import Eyegore


__all__ = ['Pipe','FlatFielding','Report','Eyegore','__version__']

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
