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

def _get_OffsetsHisto_dict(readmode,temperature):
    """ """
    OffsetsHisto_dict = dict(
        figname='BIAS_OFFSETS_HISTOS_%s_%s.png' % (readmode, temperature),
        caption='Offsets Distributions across all quadrants of the FPA.',
        meta=dict(
        title='Offsets Distribution',
        xlabel='ADU',
        ylabel='N',
        doLegend=True,
        corekwargs=dict(pre=dict(linestyle='-', color='r'),
            img=dict(linestyle='--', color='g'),
            ove=dict(linestyle='-.', color='b'))))

    return OffsetsHisto_dict

def _get_DiffOffsetsMap_dict(readmode, temperature):
    """ """

    DiffOffsetsMap_dict = dict(
        figname='BIAS_%s_%s_DiffOffsetsMap.png' % (readmode, temperature),
        caption='Difference in Offset with Reference image [ADU].',
        meta=dict(suptitle='Difference in Offsets',
                    ColorbarText='ADU'))

    return DiffOffsetsMap_dict

def _get_RatioRonMap_dict(readmode, temperature):
    """ """

    RatioRonMap_dict = dict(
        figname='BIAS_%s_%s_RatioRonMap.png' % (readmode, temperature),
        caption='Ratio of RON measured / reference.',
        meta=dict(suptitle='Ratio RONs',
                    ColorbarText='Ratio Measured/Reference'))

    return RatioRonMap_dict


def _get_VPROFS_dict(readmode,temperature):
    """ """

    VPROFS_dict = dict(
    figname='BIAS_%s_%s_VPROFS.png' % (readmode, temperature),
    caption='Vertical Average Profiles (along columns). Median value across pre/img/ove profiles has been subtracted.',
    meta=dict(
        title='Avg. Column Offset Profiles.',
        doLegend=True,
        ylim=[-20,20],
        ylabel=r'$\Delta ADU$',
        xlabel='row nr.',
        corekwargs=dict(pre=dict(marker='.', linestyle='', color='r', ms=3),
            img=dict(marker='.', linestyle='', color='g', ms=3),
            ove=dict(marker='.', linestyle='', color='b', ms=3))))

    return VPROFS_dict

def _get_HPROFS_dict(readmode,temperature):
    """ """

    HPROFS_dict = dict(
    figname='BIAS_%s_%s_HPROFS.png' % (readmode, temperature),
    caption='Horizontal Average Profiles (along rows). Split in pre-scan-to-image (left, red) and'+\
            ' image-to-over-scan cuts (right, blue). Image area boundaries indicated by vertical lines.',
    meta=dict(
        title='Avg. Row Offset Profiles.',
        doLegend=True,
        ylim=[-10,20],
        ylabel=r'$\Delta ADU$',
        xlabel='pixel',
        corekwargs=dict(pre=dict(marker='.', linestyle='', color='r', ms=3),
            pre2img=dict(marker='', linestyle='--', color='r'),
            ove=dict(marker='.', linestyle='', color='b', ms=3),
            img2ove=dict(marker='', linestyle='--', color='b'),
            )))

    return HPROFS_dict



def get_Bfigs(readmode, temperature):
    """ """

    Bfigs = dict()
    Bfigs['raw_img'] = [figclasses.Fig_Dynamic, \
            _get_raw_img_dict(readmode, temperature)]
    Bfigs['OFFSETS_HISTO'] = [figclasses.Fig_Dynamic, 
            _get_OffsetsHisto_dict(readmode,temperature)]
    Bfigs['DIFFOFFSETSMAP'] = [figclasses.Fig_Dynamic, 
            _get_DiffOffsetsMap_dict(readmode,temperature)]
    Bfigs['RATIORONSMAP'] = [figclasses.Fig_Dynamic, 
            _get_RatioRonMap_dict(readmode,temperature)]
    Bfigs['VPROFS'] = [figclasses.Fig_Dynamic, 
            _get_VPROFS_dict(readmode,temperature)]
    Bfigs['HPROFS'] = [figclasses.Fig_Dynamic, 
            _get_HPROFS_dict(readmode,temperature)]
    Bfigs['BlueScreen'] = [figclasses.BlueScreen, dict()]

    return Bfigs


def get_CDP_lib(readmode,temperature):
    """ """

    OFF_cdp = cdp.Json_CDP()
    OFF_cdp.rootname = 'BIAS_%s_%s_OFFSETS' % (readmode, temperature)
    RON_cdp = cdp.Json_CDP()
    RON_cdp.rootname = 'BIAS_%s_%s_RONS' % (readmode, temperature)

    CDP_lib = dict(OFF_JSON=OFF_cdp,
                    RON_JSON=RON_cdp)

    return CDP_lib