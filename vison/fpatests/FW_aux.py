"""

Auxiliary Functions and resources to FWD_WARM.

Created on Tue Oct 29 14:20:00 2019

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

FW_img_dict = dict(
        figname='FW_Img.png',
        caption='PENDING CAPTION')

FW_Ramps_dict = dict(
        figname='FW_Ramps.png',
        caption='PENDING CAPTION',
        meta=dict(doLegend=True,
                  corekwargs=dict(E=dict(marker='',linestyle='-',color='r'),
                                  F=dict(marker='',linestyle='-',color='g'),
                                  G=dict(marker='',linestyle='-',color='b'),
                                  H=dict(marker='',linestyle='-',color='m'))))


def get_FWfigs():
    """ """
    
    FWfigs = dict()
    
    FWfigs['FW_img'] = [figclasses.Fig_Dynamic, FW_img_dict]
    FWfigs['FW_RAMPS'] = [figclasses.Fig_Dynamic, FW_Ramps_dict]
    FWfigs['BlueScreen'] = [figclasses.BlueScreen, dict()]
    
    return FWfigs


def get_CDP_lib():
    """ """
    
    HER_cdp = cdp.Tables_CDP()
    HER_cdp.rootname = 'HER_TB'
    
    RSLOPE_cdp = cdp.Tables_CDP()
    RSLOPE_cdp.rootname='RSLOPE_TB'
    
    CDP_lib = dict(HER_TB=HER_cdp,
                   RSLOPE_TB=RSLOPE_cdp)
    
    return CDP_lib