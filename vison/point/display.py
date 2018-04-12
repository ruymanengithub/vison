# -*- coding: utf-8 -*-
"""

Display Library for Point-Source Analysis
=========================================

:requires: matplotlib
:author: Ruyman Azzollini
:contact: r.azzollini _at_ ucl.ac.uk

Created on Fri Apr 21 14:02:57 2017

"""
# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
from vison.pipe import lib as pilib
from vison.point import lib as polib

import matplotlib
matplotlib.rc('text', usetex=True)
matplotlib.rcParams['font.size'] = 17
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('axes', linewidth=1.1)
matplotlib.rcParams['legend.fontsize'] = 11
matplotlib.rcParams['legend.handlelength'] = 3
matplotlib.rcParams['xtick.major.size'] = 5
matplotlib.rcParams['ytick.major.size'] = 5
import matplotlib.pyplot as plt

# END IMPORT


def show_spots_allCCDs(spots_bag, title='', filename='', dobar=True):
    """ """

    Quads = ['E', 'F', 'H', 'G']  # not a typo

    fkpointnames = [['ALPHA', 'NONE', 'BRAVO'],
                    ['NONE', 'CHARLIE', 'NONE'],
                    ['DELTA', 'NONE', 'ECHO']]

    nx, ny = spots_bag[spots_bag.keys()[0]][Quads[0]]['ALPHA']['stamp'].shape

    NX = nx*3*2
    NY = ny*3*2

    imgs = []

    for CCDix in range(1, 4):

        iimg = np.zeros((NX, NY), dtype='float32')

        for i in range(2):
            for j in range(2):
                Q = Quads[i*2+j]
                for k in range(3):
                    for l in range(3):
                        spotname = fkpointnames[k][l]
                        if spotname != 'NONE':
                            stamp = spots_bag['CCD%i' %
                                              CCDix][Q][spotname]['stamp'].copy()
                            x0 = nx*(i*3+k)
                            x1 = x0+nx
                            y0 = ny*(j*3+l)
                            y1 = y0+ny
                            iimg[x0:x1, y0:y1] = stamp.copy()
        imgs.append(iimg)

    fig = plt.figure(figsize=(12, 6))
    axs = []

    for CCDix in range(1, 4):

        axs.append(fig.add_subplot(1, 3, CCDix))
        axs[-1].imshow(imgs[CCDix-1])
        axs[-1].set_title('CCD%i' % CCDix)

    plt.tight_layout()
    plt.suptitle(title)

    if filename == '':
        plt.show()
    else:
        plt.savefig(filename)
