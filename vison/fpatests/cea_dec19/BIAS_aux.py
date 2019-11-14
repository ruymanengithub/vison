"""

Auxiliary Functions and resources to FPA_BIAS.

Created on Thu Nov 14 10:50:00 2019

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

def _get_raw_img_dict(readmode, temperature):

    raw_img_dict = dict(
    figname='BIAS_Raw_Img_%s_%s.png' % (readmode, temperature),
    caption='FPA Bias Image Display: %s - %s' % (readmode, temperature),
    meta=dict(suptitle='%s - %s' % (readmode, temperature)))

    return raw_img_dict



def get_Bfigs(readmode, temperature):
    """ """

    Bfigs = dict()
    Bfigs['raw_img'] = [figclasses.Fig_Dynamic, _get_raw_img_dict(readmode, temperature)]
    Bfigs['BlueScreen'] = [figclasses.BlueScreen, dict()]

    return Bfigs


def get_CDP_lib():
    """ """


    CDP_lib = dict()

    return CDP_lib