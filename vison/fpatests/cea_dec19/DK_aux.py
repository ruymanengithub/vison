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

DK_img_dict = dict(
    figname='DARK_Img.png',
    caption='FPA Image: "Dark" (offset subtracted).',
    meta=dict(suptitle='FPA Image'))

DK_Histo_dict = dict(
        figname='DARK_HISTOS.png',
        caption='Distribution of Dark Current estimates across all quadrants of the FPA [e-/pix/hr]',
        meta=dict(
            title='Dark Current Distribution',
            xlabel='[e-/pix/hr]',
            ylabel='N',
            doLegend=False,
            corekwargs=dict(linestyle='-', color='r')))


DK_VPROFS_dict = dict(
    figname='DARK_VPROFS.png',
    caption='Vertical (along columns) Average Profiles. Median value across pre/img/ove profiles subtracted.',
    meta=dict(
        title='Avg. Column Profiles',
        doLegend=True,
        corekwargs=dict(pre=dict(marker='.', linestyle='', color='r', ms=3),
            img=dict(marker='.', linestyle='', color='g', ms=3),
            ove=dict(marker='.', linestyle='', color='b', ms=3))))

DK_HPROFS_dict = dict(
    figname='DARK_HPROFS.png',
    caption='Horizontal Average Profiles (along rows). Split in pre-scan-to-image (left, red) and'+\
            ' image-to-over-scan cuts (right, blue). Image area boundaries indicated by vertical lines.',
    meta=dict(
        title='Avg. Row Profiles.',
        doLegend=True,
        ylim=[-10,20],
        corekwargs=dict(pre=dict(marker='.', linestyle='', color='r', ms=3),
            pre2img=dict(marker='', linestyle='--', color='r'),
            ove=dict(marker='.', linestyle='', color='b', ms=3),
            img2ove=dict(marker='', linestyle='--', color='b'),
            )))

DKMap_dict = dict(
        figname='DARK_Map.png',
        caption='Dark Current rates across FPA in e-/pix/hour.',
        meta=dict(suptitle='Dark Current Rates',
                    ColorbarText='e-/pix/hour'))

def get_DKfigs():
    """ """

    DKfigs = dict()
    DKfigs['DK_img'] = [figclasses.Fig_Dynamic, DK_img_dict]
    DKfigs['DK_HISTO'] = [figclasses.Fig_Dynamic,DK_Histo_dict]
    DKfigs['VPROFS'] = [figclasses.Fig_Dynamic,DK_VPROFS_dict]
    DKfigs['HPROFS'] = [figclasses.Fig_Dynamic,DK_HPROFS_dict]
    DKfigs['DKMAP'] = [figclasses.Fig_Dynamic,DKMap_dict]
    DKfigs['BlueScreen'] = [figclasses.BlueScreen, dict()]

    return DKfigs


def get_CDP_lib():
    """ """

    CDP_lib = dict()

    return CDP_lib