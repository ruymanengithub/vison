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
    caption='PENDING CAPTION',
    meta=dict(suptitle='FPA Image'))


def get_CIfigs():
    """ """

    CIfigs = dict()

    CIfigs['CI_img'] = [figclasses.Fig_Dynamic, CI_img_dict]
    CIfigs['BlueScreen'] = [figclasses.BlueScreen, dict()]

    return CIfigs


def get_CDP_lib():
    """ """

    CDP_lib = dict()

    return CDP_lib