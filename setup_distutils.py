"""

:History:

22 Dec 2016: created

:author: Ruyman Azzollini (MSSL)
:contact: r.azzollini_at_ucl.ac.uk
"""

from setuptools import setup, Extension
#from numpy.distutils.core import setup,Extension

import sys

from pdb import set_trace as stop



setup(
    name="vison",
    version="0.1",
    description="Euclid VIS Ground Calibration Pipeline",
    author="Ruyman Azzollini",
    author_email="r.azzollini_at_ucl.ac.uk",
    url="https://github.com/ruymanengithub/vison",
    long_description=__doc__,
    packages=['vison',
       'vison.pipeline',
       'vison.flat',
       'vison.point',
       'vison.inject',
       'vison.pump',
       'vison.datamodel',
       'vison.sandbox',
       'vison.support'],
    package_dir={'vison':'vison/',
    'vison.pipeline':'vison/pipeline/',
    'vison.flat':'vison/flat/',
    'vison.point':'vison/point/',
    'vison.inject':'vison/inject/',
    'vison.pump':'vison/pump/',
    'vison.datamodel':'vison/datamodel/',
    'vison.sandbox':'vison/sandbox/',
    'vison.support':'vison/support/'},
    package_data={'vison':['vison/data/*']},
#    scripts=['vison/scripts/vison'],
    include_package_data=True,
    zip_safe=False,
)
