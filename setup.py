"""

:History:

22 Dec 2016: created

:author: Ruyman Azzollini (MSSL)
:contact: r.azzollini_at_ucl.ac.uk

"""

from numpy.distutils.core import setup,Extension
from numpy.distutils.misc_util import Configuration

import sys
import os

from pdb import set_trace as stop


def configuration():
    
    config = Configuration()
    
    config.add_subpackage('vison')
    config.add_subpackage('vison/analysis')
    config.add_subpackage('vison/datamodel')
    config.add_subpackage('vison/eyegore')
    config.add_subpackage('vison/flat')
    config.add_subpackage('vison/inject')
    config.add_subpackage('vison/pipe')
    config.add_subpackage('vison/point')
    config.add_subpackage('vison/pump')
    config.add_subpackage('vison/sandbox')
    config.add_subpackage('vison/support')
    
    #config.add_subpackage('vison/doc')
        
    config.add_data_dir(['vison/data','vison/data'])
    config.add_data_dir(['vison/doc','vison/doc'])  
    
    config.add_scripts(['vison/scripts/HKmonitor.py',
           'vison/scripts/QLA.py',
           'vison/scripts/quickds9.py'])
   
    
    return config


def setup_package():
        
    metadata = dict(
     name="vison",
     version="0.1",
     description="Euclid VIS Ground Calibration Pipeline",
     author="Ruyman Azzollini",
     author_email="r.azzollini_at_ucl.ac.uk",
     url="https://github.com/ruymanengithub/vison",
     configuration=configuration)
    
    setup(**metadata)
    
    return
    

if __name__ == '__main__': setup_package()

