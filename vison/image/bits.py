#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

NEEDSREVISION

Image bits analysis tools.

Created on Thu Sep 14 15:54:14 2017

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np

import os
import string as st

from vison.datamodel import ccd as ccdmod

import matplotlib
matplotlib.rc('text', usetex=True)
matplotlib.rcParams['font.size'] = 17
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('axes', linewidth=1.1)
matplotlib.rcParams['legend.fontsize'] = 11
matplotlib.rcParams['legend.handlelength'] = 3
matplotlib.rcParams['xtick.major.size'] = 5
matplotlib.rcParams['ytick.major.size'] = 5
from matplotlib import cm
import matplotlib.pyplot as plt

# END IMPORT

isthere = os.path.exists

Quads = ccdmod.Quads


def show_histo_adus(qdata, title=''):

    bins = np.arange(0, 2**16-1)
    histo = np.histogram(qdata, bins=bins, density=True)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    PDF = histo[0]
    bins = histo[1]

    bincenters = (bins[0:-1]+bins[1:])/2.
    #binwidth = bins[1]-bins[0]

    ax.semilogx(bincenters, PDF, color='blue')

    ax.set_title(title)

    ax.set_xlabel('ADU')
    ax.set_ylabel('Freq.')
    ax.set_xlim([1, 2**16])

    plt.tight_layout()
    plt.show()
    plt.close()


# def extract_bits(integer):
#    #bits = [int(item) for item in list("{0:b}".format(int(integer)).rjust(16,'0'))]
#    bits = np.array(list(bin(int(integer))[2:].rjust(16,'0'))).astype('int')[::-1]
#    return bits


def show_histo_bits(ccdobj, suptitle='', figname=''):

    pQuads = ['E', 'F', 'H', 'G']

    bins = np.arange(0, 16)+0.5

    axs = []
    fig = plt.figure()

    for iQ, pQ in enumerate(pQuads):

        axs.append(fig.add_subplot(2, 2, iQ+1))
        qdata = ccdobj.get_quad(pQ, canonical=True)

        vector = qdata.astype('int32').flatten()
        bits = (((vector[:, None] & (1 << np.arange(16)))) > 0).astype(int)
        bitsmean = bits.mean(axis=0)

        axs[-1].bar(bins, bitsmean, width=1,
                    edgecolor='blue', facecolor='white')
        axs[-1].axhline(y=0.5, ls='--', color='r')

        axs[-1].set_xticks(bins)
        axs[-1].set_xticklabels(np.arange(16).astype('S'))

        # ax.plot(bins,bitsmean,color='blue')

        axs[-1].set_title(pQ)

        axs[-1].set_xlabel('bit')
        axs[-1].set_ylabel('Average')

    plt.suptitle(suptitle)
    plt.tight_layout()
    plt.subplots_adjust(top=0.85)

    if figname != '':
        plt.savefig(figname)
    else:
        plt.show()
    plt.close()


def _get_corrmap(bits):
    """ """

    Nbits = bits.shape[1]
    corrmx = np.zeros((Nbits, Nbits), dtype='float32')

    mbits = (bits - bits.mean(axis=0))
    sbits = np.std(mbits, axis=0)

    for i in range(Nbits):
        for j in range(i, Nbits):

            EXY = (mbits[:, i]*mbits[:, j]).mean()
            sx = sbits[i]
            sy = sbits[j]
            corrmx[i, j] = EXY/(sx*sy)

    corrmx[np.isnan(corrmx)] = 0.
    corrmx = corrmx + corrmx.T - np.diag(np.diag(corrmx))

    return corrmx


def show_correlation_bits(ccdobj, suptitle='', figname=''):

    pQuads = ['E', 'F', 'H', 'G']

    bins = np.arange(0, 16, 2)+0.5

    axs = []
    fig = plt.figure()

    for iQ, pQ in enumerate(pQuads):

        axs.append(fig.add_subplot(2, 2, iQ+1))
        qdata = ccdobj.get_quad(pQ, canonical=True)

        vector = qdata.astype('int32').flatten()
        bits = (((vector[:, None] & (1 << np.arange(16)))) > 0).astype('int8')
        corrmx = _get_corrmap(bits)

        #xx,yy = np.meshgrid(np.arange(16),np.arange(16))
        #corrmx = 1.*xx/xx.max()
        #corrmx -= corrmx.mean()

        im = axs[-1].imshow(corrmx, origin='lower left', cmap=cm.rainbow)

        axs[-1].set_xticks(bins)
        axs[-1].set_xticklabels(np.arange(0, 16, 2).astype('S'))

        axs[-1].set_yticks(bins)
        axs[-1].set_yticklabels(np.arange(0, 16, 2).astype('S'))

        axs[-1].set_title(pQ)

    cax = fig.add_axes([0.85, 0.1, 0.03, 0.8])
    fig.colorbar(im, cax=cax)

    plt.suptitle(suptitle)
    plt.tight_layout()
    plt.subplots_adjust(top=0.85)
    plt.subplots_adjust(right=0.85)

    if figname != '':
        plt.savefig(figname)
    else:
        plt.show()
    plt.close()


def test():

    # data = dict(fits='XTALK_July2017/18July2017/Closed_ROE_BOX/With_Lid/MODs/Inductor_shields/L523_AND_L516_AND_L517_shielded/Channel_3/H/EUC_6_180717d170933T_ROE1_CCD2.fits',
    #            Q='F')
    data = dict(fits='PIXBOUNCE_images_atDisk3/17_7_17/FPGA_568C/FWD/Baseline/Channel_1/17_07/EUC_2_170717D093830T_ROE1_CCD1.fits', Q='F')

    # data = dict(fits = 'PIXBOUNCE_images_atDisk3/17_7_17/FPGA_568C/RWD/Channel1/17_07/EUC_1_170717D162109T_ROE1_CCD1.fits',
    #            Q='E') # positive

    fits = data['fits']
    ccdobj = ccdmod.CCD(fits)
    # Q=data['Q']

    suptitle = '%s' % (st.replace(os.path.split(fits)[-1], '_', '\_'))
    # show_histo_adus(qdata,title=title)

    show_histo_bits(ccdobj, suptitle=suptitle)
