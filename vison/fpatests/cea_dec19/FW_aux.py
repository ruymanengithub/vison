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
    caption='FPA FWD WARM Image Display',
    meta=dict(suptitle='FPA Image - FWD-WARM'))

FW_Ramps_dict = dict(
    figname='FW_Ramps.png',
    caption='Avg. Column Fluence Profiles [ADU].',
    meta=dict(
        suptitle='Avg. Column Fluence Profiles.',
        doLegend=True,
        corekwargs=dict(E=dict(marker='', linestyle='-', color='r'),
            F=dict(marker='', linestyle='-', color='g'),
            G=dict(marker='', linestyle='-', color='b'),
            H=dict(marker='', linestyle='-', color='m'))))

FW_SlopesMap_dict = dict(
    figname='FW_SlopesMap.png',
    caption='Ramp Slopes in ADU per pixel row.',
    meta=dict(suptitle='Ramp Slopes',
                ColorbarText='ADU/row'))

FW_DiffSlopesMap_dict = dict(
    figname='FW_DiffSlopesMap.png',
    caption='Difference in ramp slope with reference image (from GRCALCAMP),'+\
            ' in ADU per pixel row.',
    meta=dict(suptitle='Diff. Ramp Slopes',
                ColorbarText=r'$\Delta$ADU/row'
            ))

FW_HERPROFS_dict = dict(
    figname='FW_HERPROFS.png',
    caption='Hard Edge Response curves in the image to over-scan transition region, '+\
            'all 144 quadrants.',
    meta=dict(
        title='H.E.R. Curves',
        xlabel='Pixel',
        ylabel='HER [frac]',
        ylim=[-2.E-4, 5.e-4],
        xlim=[9, 15],
        corekwargs=dict(linestyle='-', marker='')))

FW_HERvalsMap_dict = dict(
    figname='FW_HERvalsMap.png',
    caption='HER values (adim.) in 1st dark pixel across FPA.',
    meta=dict(suptitle='HER (1st pixel) Map',
                ColorbarText='Fraction'))


FW_BITHISTOS_dict = dict(
    figname='FW_BITS_Histos.png',
    caption='Bit Population histograms across the 144 FPA quadrants (different colours).',
    meta=dict(
        title='Bit Histograms',
        xlabel='Bit',
        ylabel='Rel. Freq.',
        ylim=[0., 1.1],
        xlim=[0, 16],
        corekwargs=dict(linestyle=' ', marker='.')))

FW_DiffOffsetsMap_dict = dict(
        figname='FW_DiffOffsetsMap.png',
        caption='Difference in Offset with Reference image [ADU].',
        meta=dict(suptitle='Difference in Offsets',
                    ColorbarText='ADU'))

FW_RatioRonMap_dict = dict(
        figname='FW_RatioRonMap.png',
        caption='Ratio of RON measured / reference.',
        meta=dict(suptitle='Ratio RONs',
                    ColorbarText='Ratio Measured/Reference'))



def get_FWfigs():
    """ """

    FWfigs = dict()

    FWfigs['FW_img'] = [figclasses.Fig_Dynamic, FW_img_dict]
    FWfigs['FW_RAMPS'] = [figclasses.Fig_Dynamic, FW_Ramps_dict]
    FWfigs['SLOPESMAP'] = [figclasses.Fig_Dynamic, FW_SlopesMap_dict]
    FWfigs['DIFFSLOPESMAP'] = [figclasses.Fig_Dynamic, FW_DiffSlopesMap_dict]
    FWfigs['DIFFOFFSETSMAP'] = [figclasses.Fig_Dynamic, FW_DiffOffsetsMap_dict]
    FWfigs['RATIORONSMAP'] = [figclasses.Fig_Dynamic, FW_RatioRonMap_dict]
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
