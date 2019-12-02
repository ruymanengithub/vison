"""

Auxiliary Functions and resources to FPA_CHINJ.

Created on Thu Nov 14 13:07:00 2019

:author: Ruyman Azzollini

"""


# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import os
from collections import OrderedDict

from vison.plot import figclasses
from vison.datamodel import cdp


# END IMPORT

CI_img_dict = dict(
    figname='CHINJ_Img.png',
    caption='FPA Image Display: Charge Injection',
    meta=dict(suptitle='FPA Image'))

INJ_Profs_dict = dict(
    figname='CHINJ_Injection_Profiles.png',
    caption='Average Injection [ADU] profiles across the FPA.',
    meta=dict(
            suptitle='Avg. Injection Profiles.',
            doLegend=True,
            corekwargs=dict(E=dict(marker='', linestyle='-', color='r'),
                            F=dict(marker='', linestyle='-', color='g'),
                            G=dict(marker='', linestyle='-', color='b'),
                            H=dict(marker='', linestyle='-', color='m'))))

INJ_Diff_dict = dict(
    figname='CHINJ_DiffInjectionMap.png',
    caption='Difference in injection levels [ADU] with reference image from GRCALCAMP.',
    meta=dict(suptitle='Diff. Injection Levels'))


def get_CIfigs():
    """ """

    CIfigs = dict()

    CIfigs['CI_img'] = [figclasses.Fig_Dynamic, CI_img_dict]
    CIfigs['INJ_PROFS'] = [figclasses.Fig_Dynamic, INJ_Profs_dict]
    CIfigs['DIFF_INJMAP'] = [figclasses.Fig_Dynamic, INJ_Diff_dict]
    CIfigs['BlueScreen'] = [figclasses.BlueScreen, dict()]

    return CIfigs


def get_CDP_lib():
    """ """

    INJ_cdp = cdp.Tables_CDP()
    INJ_cdp.rootname = 'INJ_TB'

    CDP_lib = dict(INJ_TB=INJ_cdp)

    return CDP_lib