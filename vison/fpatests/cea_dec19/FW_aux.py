"""

Auxiliary Functions and resources to FWD_WARM.

Created on Tue Oct 29 14:20:00 2019

:author: Ruyman Azzollini

"""


# IMPORT STUFF
# from pdb import set_trace as stop
# import numpy as np
# import os
# from collections import OrderedDict

from vison.plot import figclasses
from vison.datamodel import cdp


# END IMPORT

FW_img_dict = dict(
    figname='FW_Img.png',
    caption='PENDING CAPTION',
    meta=dict(suptitle='FPA Image'))

FW_Ramps_dict = dict(
    figname='FW_Ramps.png',
    caption='PENDING CAPTION',
    meta=dict(
        suptitle='Avg. Column Fluence Profiles.',
        doLegend=True,
        corekwargs=dict(E=dict(marker='', linestyle='-', color='r'),
            F=dict(marker='', linestyle='-', color='g'),
            G=dict(marker='', linestyle='-', color='b'),
            H=dict(marker='', linestyle='-', color='m'))))

FW_SlopesMap_dict = dict(
    figname='FW_SlopesMap.png',
    caption='PENDING CAPTION',
    meta=dict(suptitle='Ramp Slopes'))

FW_DiffSlopesMap_dict = dict(
    figname='FW_DiffSlopesMap.png',
    caption='PENDING CAPTION',
    meta=dict(suptitle='Diff. Ramp Slopes'))

FW_HERPROFS_dict = dict(
    figname='FW_HERPROFS.png',
    caption='PENDING CAPTION',
    meta=dict(
        title='H.E.R. Curves',
        xlabel='Pixel',
        ylabel='HER [frac]',
        ylim=[-2.E-4, 5.e-4],
        xlim=[9, 15],
        corekwargs=dict(linestyle='-', marker='')))

FW_HERvalsMap_dict = dict(
    figname='FW_HERvalsMap.png',
    caption='PENDING CAPTION',
    meta=dict(suptitle='HER (1st pixel) Map'))


FW_BITHISTOS_dict = dict(
    figname='FW_BITS_Histos.png',
    caption='PENDING CAPTION',
    meta=dict(
        title='Bit Histograms',
        xlabel='Bit',
        ylabel='Rel. Freq.',
        ylim=[0., 1.1],
        xlim=[0, 16],
        corekwargs=dict(linestyle=' ', marker='.')))


def get_FWfigs():
    """ """

    FWfigs = dict()

    FWfigs['FW_img'] = [figclasses.Fig_Dynamic, FW_img_dict]
    FWfigs['FW_RAMPS'] = [figclasses.Fig_Dynamic, FW_Ramps_dict]
    FWfigs['SLOPESMAP'] = [figclasses.Fig_Dynamic, FW_SlopesMap_dict]
    FWfigs['DIFFSLOPESMAP'] = [figclasses.Fig_Dynamic, FW_DiffSlopesMap_dict]
    FWfigs['HERPROFS'] = [figclasses.Fig_Dynamic, FW_HERPROFS_dict]
    FWfigs['HERVALSMAP'] = [figclasses.Fig_Dynamic, FW_HERvalsMap_dict]
    FWfigs['BITHISTOS'] = [figclasses.Fig_Dynamic, FW_BITHISTOS_dict]
    FWfigs['BlueScreen'] = [figclasses.BlueScreen, dict()]

    return FWfigs


def get_CDP_lib():
    """ """

    HER_cdp = cdp.Tables_CDP()
    HER_cdp.rootname = 'HER_TB'

    RSLOPE_cdp = cdp.Tables_CDP()
    RSLOPE_cdp.rootname = 'RSLOPE_TB'

    CDP_lib = dict(HER_TB=HER_cdp,
                   RSLOPE_TB=RSLOPE_cdp)

    return CDP_lib
