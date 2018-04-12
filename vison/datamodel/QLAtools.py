# -*- coding: utf-8 -*-

"""

Quick-Look-Analysis Tools.

:History:
Created on Wed Mar 16 11:31:58 2016

@author: Ruyman Azzollini
"""

# IMPORT STUFF
import numpy as np
import os
from matplotlib import pyplot as plt
import string as st

from vison.datamodel import ccd
from vison.support import latex

from pdb import set_trace as stop
# END IMPORT

isthere = os.path.exists

QUADS = ['E', 'F', 'G', 'H']


latex_templates = {'table_stats':
                   ['\\begin{center}',
                    '\\resizebox{0.3\\textwidth}{!}{',
                    '\\begin{tabular}{|c|c|c|c|c|c|c|c|}',
                    '\hline',
                    '\multicolumn{7}{c}{%s}\\\\'
                    '\hline',
                    'section & mean & std & min & max & p25 & p50 & p75\\\\',
                    '\hline',
                    'prescan & %.2f & %.3f & %.2f & %.2f & %.2f & %.2f & %.2f\\\\',
                    'image & %.2f & %.3f & %.2f & %.2f & %.2f & %.2f & %.2f\\\\',
                    'overscan & %.2f & %.3f & %.2f & %.2f & %.2f & %.2f & %.2f\\\\',
                    '\hline',
                    '\end{tabular}}',
                    '\end{center}']}


def getacrosscolscut(CCDobj):
    """ """

    result = {}

    for iQ, QUAD in enumerate(QUADS):

        data = CCDobj.get_quad(QUAD).copy()

        ncols = data.shape[0]

        stats = np.zeros((2, ncols))

        for icol in range(ncols):
            col = data[icol, :]
            stats[:, icol] = (np.mean(col), np.std(col))

        result[QUAD] = {}
        result[QUAD] = stats.copy()

    return result


def getacrossrowscut(CCDobj):
    """ """

    result = {}

    for iQ, QUAD in enumerate(QUADS):

        prestart, preend, imgstart, imgend, ovstart, ovend = CCDobj.getsectioncollims(
            QUAD)

        prestart += 3
        preend -= 3
        imgstart += 3
        imgend -= 3
        ovstart += 3
        ovend -= 3

        data = CCDobj.get_quad(QUAD).copy()

        nrows = data.shape[1]

        prestats = np.zeros((2, nrows))
        ovstats = np.zeros((2, nrows))
        imgstats = np.zeros((2, nrows))

        for irow in range(nrows):

            prescan = data[prestart:preend, irow]
            ovscan = data[ovstart:ovend, irow]
            imgline = data[imgstart:imgend, irow]

            prestats[:, irow] = (np.mean(prescan), np.std(prescan))
            ovstats[:, irow] = (np.mean(ovscan), np.std(ovscan))
            imgstats[:, irow] = (np.mean(imgline), np.std(imgline))

        result[QUAD] = {}
        result[QUAD]['prescan'] = prestats.copy()
        result[QUAD]['overscan'] = ovstats.copy()
        result[QUAD]['image'] = imgstats.copy()

    return result


def getsectionstats(CCDobj, QUAD, section, xbuffer=(0, 0), ybuffer=(0, 0)):
    """ """

    prestart, preend, imgstart, imgend, ovstart, ovend = CCDobj.getsectioncollims(
        QUAD)

    if section == 'prescan':
        x0, x1 = prestart, preend-xbuffer[1]
    elif section == 'overscan':
        x0, x1 = ovstart, ovend
    elif section == 'image':
        x0, x1 = imgstart, imgend

    x0 += xbuffer[0]
    x1 -= xbuffer[1]

    quaddata = CCDobj.get_quad(QUAD)
    NX, NY = quaddata.shape

    y0, y1 = (0, NY-1)
    y0 += ybuffer[0]
    y1 -= ybuffer[1]

    subquad = quaddata[x0:x1, y0:y1]

    stats = [np.mean(subquad), np.std(subquad), np.min(subquad), np.max(subquad),
             np.percentile(subquad, q=25), np.percentile(subquad, q=50), np.percentile(subquad, q=75)]

    return stats


def plotQuads(CCDobj, filename=None, suptitle=''):
    """ """
    figP = plt.figure(figsize=(6, 6))

    figP.suptitle(suptitle)
    axP = []

    pQUADS = ['E', 'F', 'H', 'G']

    for iQ, QUAD in enumerate(pQUADS):

        data = CCDobj.get_quad(QUAD).copy()
        ixnozero = np.where(data != 0.)
        if len(ixnozero[0] != 0):
            minP = np.percentile(data[ixnozero], q=30)
            maxP = np.percentile(data[ixnozero], q=70)
        else:
            minP = 0
            maxP = 0

        axP.append(figP.add_subplot(2, 2, iQ+1))
        axP[-1].imshow(data.transpose(), origin='lower left', cmap='hot',
                       clim=(minP, maxP))
        axP[-1].set_title(QUAD)

    plt.tight_layout()

    if filename is not None:
        plt.savefig(filename)
    else:
        plt.show()
    plt.close('all')


def plotAcROWcuts(dissection, filename=None, suptitle=''):
    """ """

    figP = plt.figure(figsize=(6, 6))
    figP.suptitle(suptitle)

    axP = []

    pQUADS = ['E', 'F', 'H', 'G']

    for iQ, QUAD in enumerate(pQUADS):

        data = dissection[QUAD].copy()
        ncols = data.shape[1]

        col = np.arange(ncols)+1
        meanval = data[0, :].copy()

        axP.append(figP.add_subplot(2, 2, iQ+1))
        axP[-1].plot(col, meanval)
        axP[-1].set_title(QUAD)

        if QUAD in ['E', 'H']:
            axP[-1].set_ylabel('ADU')
        if QUAD in ['H', 'G']:
            axP[-1].set_xlabel('col')

        axP[-1].ticklabel_format(style='sci', axis='x', scilimits=(2, 2))

    plt.tight_layout()

    if filename is not None:
        plt.savefig(filename)
    else:
        plt.show()
    plt.close('all')


def plotAcCOLcuts(dissection, filename=None, suptitle=''):
    """ """

    figP = plt.figure(figsize=(6, 6))
    figP.suptitle(suptitle)

    axP = []

    pQUADS = ['E', 'F', 'H', 'G']

    colors = ['b', 'g', 'r']

    for iQ, QUAD in enumerate(pQUADS):

        data = dissection[QUAD]

        axP.append(figP.add_subplot(2, 2, iQ+1))

        for ia, area in enumerate(['prescan', 'image', 'overscan']):

            nrows = data[area].shape[1]
            row = np.arange(nrows)+1
            meanval = data[area][0, :].copy()

            axP[-1].plot(row, meanval, color=colors[ia])

        axP[-1].set_title(QUAD)

        if QUAD in ['E', 'H']:
            axP[-1].set_ylabel('ADU')
        if QUAD in ['H', 'G']:
            axP[-1].set_xlabel('row')

        axP[-1].ticklabel_format(style='sci', axis='x', scilimits=(2, 2))

    plt.tight_layout()

    if filename is not None:
        plt.savefig(filename)
    else:
        plt.show()
    plt.close('all')


def dissectFITS(FITSfile, path=''):
    """ """

    diss = {}

    CCDobj = ccd.CCD(FITSfile)

    # Header selection

    sel_header_keys = ['EXTNAME', 'PROGRAM', 'OBJECT', 'OBSID', 'DATE', 'EXPTIME',
                       'CHRG_INJ', 'TRAPPUMP', 'VSTART', 'VEND', 'WAVELENG',
                       'CCD_SN', 'ROE_SN', 'FPGA_VER', 'TOI_FLU', 'TOI_PUMP',
                       'TOI_READ', 'INVFLSHP', 'INVFLUSH', 'FLUSHES', 'R1CCD2TT',
                       'R1CCD2TB', 'IDL_V', 'IDH_V', 'IG1_V', 'IG2_V',
                       'ODCCD2_V', 'RDCCD2_V']

    diss['subheader'] = []

    for key in sel_header_keys:
        diss['subheader'].append('%s=%s' % (key, CCDobj.header[key]))

    OBSID = CCDobj.header['OBSID']
    DATE = CCDobj.header['DATE']
    rootfname = '%s_%s' % (OBSID, DATE)

    xbuffer = (3, 3)
    ybuffer = (3, 3)

    # STATISTICS : mean, std, min, max, p25, p50, p75

    diss['imgstats'] = {}
    diss['prestats'] = {}
    diss['overstats'] = {}

    for iQ, QUAD in enumerate(QUADS):

        diss['imgstats'][QUAD] = getsectionstats(
            CCDobj, QUAD, 'image', xbuffer=xbuffer, ybuffer=ybuffer)
        diss['prestats'][QUAD] = getsectionstats(
            CCDobj, QUAD, 'prescan', xbuffer=xbuffer, ybuffer=ybuffer)
        diss['overstats'][QUAD] = getsectionstats(
            CCDobj, QUAD, 'overscan', xbuffer=xbuffer, ybuffer=ybuffer)

    # Across Columns profiles of mean and std

    diss['acrosscolscut'] = getacrosscolscut(CCDobj)

    accolfname = '%s_ACCOL.eps' % (rootfname,)
    accolfname = os.path.join(path, accolfname)
    plotAcROWcuts(diss['acrosscolscut'], filename=accolfname, suptitle=OBSID)

    # Across Rows profiles of mean and std

    diss['acrossrowscut'] = getacrossrowscut(CCDobj)

    acrowfname = '%s_ACROW.eps' % (rootfname,)
    acrowfname = os.path.join(path, acrowfname)
    plotAcCOLcuts(diss['acrossrowscut'], filename=acrowfname, suptitle=OBSID)

    # Quadrant Stamps

    quadsfname = '%s_STAMP.eps' % rootfname
    quadsfname = os.path.join(path, quadsfname)

    plotQuads(CCDobj, filename=quadsfname, suptitle=OBSID)
    diss['STAMP'] = quadsfname
    diss['COLS'] = accolfname
    diss['ROWS'] = acrowfname

    return diss


def reportFITS(FITSfile, outpath=''):
    """ """

    assert isthere(FITSfile)

    bareimg = os.path.split(FITSfile)[-1]
    PDFroot = '%s' % os.path.splitext(bareimg)[0]

    dissection = dissectFITS(FITSfile, path=outpath)

    report = latex.LaTeX(fontsize=10)
    figlist = []
    niceFITSname = st.replace(bareimg, '_', '\\_')
    report.body.append('%s\\\\' % niceFITSname)

    report.body.append('\n')

    report.body += ['\\begin{multicols}{2}']

    report.body.append('\\begin{verbatim}')
    for line in dissection['subheader']:
        #niceline = st.replace(line,'_','\_')
        #report.body.append('%s\\\\' % niceline)
        report.body.append(line)
    report.body.append('\\end{verbatim}')
    # report.body.append('\\newline')

    for QUAD in QUADS:

        prestats = dissection['prestats'][QUAD]
        imgstats = dissection['imgstats'][QUAD]
        overstats = dissection['overstats'][QUAD]

        alignedstats = tuple([QUAD]+prestats+imgstats+overstats)

        template_stats = st.join(latex_templates['table_stats'], '__^__')
        table_stats = template_stats % alignedstats
        table_stats = st.split(table_stats, '__^__')

        report.body += table_stats

        # getacrosscolscut(CCDobj,QUAD)
        # getacrossrowscut(CCDobj,QUAD)

    report.body.append('\columnbreak')

    report.addfigtobody(dissection['STAMP'], imgsize=7)
    report.addfigtobody(dissection['COLS'], imgsize=7)
    report.addfigtobody(dissection['ROWS'], imgsize=7)

    figlist += [dissection['STAMP'], dissection['COLS'], dissection['ROWS']]

    report.body.append('\end{multicols}')

    report.Write('%s.tex' % PDFroot)

    report.Compile2PDF('%s.tex' % PDFroot, cleanafter=True, figures=figlist)
    os.system('mv %s.pdf %s/' % (PDFroot, outpath))

    PDF = os.path.join(outpath, '%s.pdf' % PDFroot)

    return PDF
