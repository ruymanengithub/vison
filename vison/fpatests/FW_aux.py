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
        caption='PENDING Caption.')

def get_FWfigs():
    """ """
    
    FWfigs = dict()
    
    FWfigs['FW_img'] = [figclasses.Fig_Husk, FW_img_dict]
    FWfigs['BlueScreen'] = [figclasses.BlueScreen, dict()]
    
    return FWfigs


def get_CDP_lib():
    """ """
    
    HER_cdp = cdp.Tables_CDP()
    HER_cdp.rootname = 'HER_TB'
    
    CDP_lib = dict(HER_TB=HER_cdp)
    
    return CDP_lib