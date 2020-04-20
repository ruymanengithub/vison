"""

:History:

22 Dec 2016: created

:author: Ruyman Azzollini (MSSL)

"""

from numpy.distutils.core import setup, Extension
from numpy.distutils.misc_util import Configuration

import sys
import os

from pdb import set_trace as stop

import versioneer


def configuration():

    config = Configuration()

    config.add_subpackage('vison')
    config.add_subpackage('vison/analysis')
    config.add_subpackage('vison/dark')
    config.add_subpackage('vison/campaign')
    config.add_subpackage('vison/datamodel')
    config.add_subpackage('vison/eyegore')
    config.add_subpackage('vison/flat')
    config.add_subpackage('vison/fpa')    
    config.add_subpackage('vison/fpatests')
    config.add_subpackage('vison/fpatests/cea_dec19')
    config.add_subpackage('vison/image')
    config.add_subpackage('vison/inject')
    config.add_subpackage('vison/metatests')
    config.add_subpackage('vison/obsolete')
    config.add_subpackage('vison/ogse')
    config.add_subpackage('vison/other')
    config.add_subpackage('vison/pipe')
    config.add_subpackage('vison/plot')
    config.add_subpackage('vison/point')
    config.add_subpackage('vison/pump')
    config.add_subpackage('vison/roe_fft')
    config.add_subpackage('vison/sandbox')
    config.add_subpackage('vison/scripts')
    config.add_subpackage('vison/support')        
    config.add_subpackage('vison/test')
    config.add_subpackage('vison/xml')
    config.add_subpackage('vison/xtalk')

    # config.add_subpackage('vison/doc')
    
    config.add_data_dir(['vison/blocks','vison/blocks'])
    config.add_data_dir(['vison/ogse_profiles','vison/ogse_profiles'])
    config.add_data_dir(['vison/data', 'vison/data'])
    config.add_data_dir(['vison/doc', 'vison/doc'])

    config.add_scripts(['vison/scripts/HKmonitor.py',
                        'vison/scripts/QLA.py',
                        'vison/scripts/quickds9.py',
                        'vison/scripts/vis_mkscripts.py',
                        'vison/scripts/vis_genDataSet.py',
                        'vison/scripts/eyegore',
                        'vison/scripts/vison_run',
                        'vison/scripts/vis_explogs_merger.py',
                        'vison/scripts/vis_run_xtalk.py',
                        'vison/scripts/vis_ROE_LinCalib.py',
                        'vison/scripts/vis_ROETAB_LinCalib.py',
                        'vison/scripts/vis_star_finder.py',
                        'vison/scripts/vis_load_DD.py',
                        'vison/scripts/vis_cosmetics_masker.py',
                        'vison/scripts/vis_clear_space.py',
                        'vison/scripts/vis_reports_merger.py',
                        'vison/scripts/vis_mksession.py',
                        'vison/scripts/vis_dnl_extract.py',
                        'vison/scripts/vison'
                        ])

    return config


def setup_package():

    metadata = dict(
        name="vison",
        description="Euclid VIS Ground Calibration Pipeline",
        author="Ruyman Azzollini",
        author_email="r.azzollini_at_ucl.ac.uk",
        url="https://github.com/ruymanengithub/vison",
        configuration=configuration,
        version=versioneer.get_version(),
        cmdclass=versioneer.get_cmdclass())

    setup(**metadata)

    return


if __name__ == '__main__':
    setup_package()
