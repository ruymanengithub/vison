

__version__ = "0.1.0"

import matplotlib
matplotlib.use("TkAgg")


from pipe import FlatFielding
from pipe.master import Pipe
from support.report import Report
from eyegore.eyegore import Eyegore


__all__ = ['Pipe','FlatFielding','Report','Eyegore','__version__']
