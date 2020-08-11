#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Aux. to task.py

Created on Fri Nov 17 19:11:53 2017

:author: raf

"""

# IMPORT STUFF
import os
from pdb import set_trace as stop
import string as st
import numpy as np
from collections import OrderedDict
import copy

from vison.datamodel.generator import generate_Explog
from vison.datamodel import HKtools
from vison.pipe import lib as pilib
from vison.datamodel import core, ccd
from vison.support import context
# END IMPORT


def check_HK(self, HKKeys, reference='command', limits='P', tag='', doReport=False, doLog=True):
    """ """

    checkers = dict(command=HKtools.check_HK_vs_command,
                    abs=HKtools.check_HK_abs)

    report_HK = checkers[reference](
        HKKeys, self.dd, limits=limits, elvis=self.elvis)

    if doReport and self.report is not None:
        nicetag = tag.replace(' ', '\ ')
        msg_HK = ['$\\bf{HK-%s}$, Parameters not passing Tests:\n' % nicetag]
        msg_HK += [''.join(['$\\bf{\\textcolor{red}{%s}}$' % (key.replace( '_', '\_'),)
                            for key in report_HK if not report_HK[key]], ', ')]
        if len(msg_HK[-1]) == 0:
            msg_HK += ['$\\bf{ALL\ PASSED}$\n']
        self.report.add_Text('\n'.join(msg_HK))
        self.report.add_Text('\\')

    if doLog and self.log is not None:
        msg_HK = ['HK-%s, Parameters not passing Tests:\n' % tag]
        msg_HK += [', '.join([key for key in report_HK if not report_HK[key]])]
        if len(msg_HK[-1]) == 0:
            msg_HK += ['ALL PASSED']
        self.log.info(msg_HK)

    return report_HK


def add_checkHK_report(self, report_HK, tag):

    msg_HK = ['$\\bf{HK-%s}$, Listing parameters not passing Tests:\n' % tag]
    msg_HK += [', '.join(['\\bf{%s}' % (key.replace('_', '\_'),)
                        for key in report_HK if not report_HK[key]])]
    if len(msg_HK[-1]) == 0:
        msg_HK += ['$\\bf{ALL\ PASSED}$\n']

    if self.report is not None:
        self.report.add_Text('\n'.join(msg_HK))
        self.report.add_Text('\\')

    return msg_HK


def create_mockexplog(self, OBSID0=1000):
    """ """

    logdefaults = {'egse_ver': self.elvis, 'con_file': 'vis_roe_config_cotsqm_273_vn.txt',
                   'fl_rdout': 0, 'ci_rdout': 0,
                   'fpga_ver': '2AC',
                   'cdpu_clk': 0, 'v_tp_mod': 123, 's_tp_mod': 31,
                   'motr': 1,
                   'R1C1_TT': -153., 'R1C2_TT': -153., 'R1C3_TT': -153.,
                   'R1C1_TB': -153., 'R1C2_TB': -153., 'R1C3_TB': -153., }

    explog = generate_Explog(self.inputs['structure'], defaults=logdefaults,
                             elvis=self.elvis, explog=None,
                             OBSID0=OBSID0, CHAMBER=self.inputs['CHAMBER'])
    return explog


def filterexposures(self, structure, explog, OBSID_lims, colorblind=False, wavedkeys=[],
                    surrogate=''):
    """Loads a list of Exposure Logs and selects exposures from test 'test'.

    The filtering takes into account an expected structure for the
    acquisition script.

    The datapath becomes another column in DataDict. This helps dealing
    with tests that run overnight and for which the input data is in several
    date-folders.


    """

    if len(OBSID_lims) == 0:
        OBSID_lims = [explog['ObsID'][0], explog['ObsID'][-1]]

    Ncols = structure['Ncols']

    if surrogate != '':
        testkey = surrogate
    else:
        testkey = self.inputs['test']
    #testkey = structure['col001']['test']

    if not colorblind:

        Filters = np.array([structure['col%03i' % i]['wave'] for i in range(1, Ncols + 1)])
        Filter = Filters[0]
        assert np.all(np.array(Filters) == Filter)

        selbool = (explog['test'] == testkey) & \
            (explog['ObsID'] >= OBSID_lims[0]) & \
            (explog['ObsID'] <= OBSID_lims[1]) & \
            (explog['wave'] == Filter)  # TESTS
    else:

        selbool = (explog['test'] == testkey) & \
            (explog['ObsID'] >= OBSID_lims[0]) & \
            (explog['ObsID'] <= OBSID_lims[1])

    explog = explog[np.where(selbool)]

    # Assess structure

    checkreport = pilib.check_test_structure(explog, structure, CCDs=[1, 2, 3],
                                             wavedkeys=wavedkeys)

    # Labeling of exposures

    explog = self.add_labels_to_explog(explog, structure)

#    explog['label'] = np.array(['NoneNoneNone']*len(explog))
#
#    frcounter = 0
#    for ic in range(1,Ncols+1):
#        _frames = structure['col%03i' % ic]['frames']
#        #print frcounter,frcounter+_frames*3
#        explog['label'][frcounter:frcounter+_frames*3] = 'col%03i' % ic
#        frcounter += _frames*3

    return explog, checkreport


def add_labels_to_explog(self, explog, structure):
    """ """
    Ncols = structure['Ncols']

    explog['label'] = np.array(['NoneNoneNone'] * len(explog))

    frcounter = 0
    for ic in range(1, Ncols + 1):
        _frames = structure['col%03i' % ic]['frames']
        #print frcounter,frcounter+_frames*3
        explog['label'][frcounter:frcounter + _frames * 3] = 'col%03i' % ic
        frcounter += _frames * 3

    return explog


def addHKPlotsMatrix(self):
    """Adds to self.report a table-figure with HK [self.HKKeys] during test."""

    figspath = self.inputs['subpaths']['figs']
    HKKeys = copy.deepcopy(self.HKKeys)

    time = self.dd.mx['time'][:, 0].copy()
    HKdata = OrderedDict()
    HKdata['time'] = time.copy()
    for key in HKKeys:
        HKdata[key] = self.dd.mx['HK_%s' % key][:].copy()

    HKfigs = []

    for key in HKKeys:

        HKlims = HKtools.HKlims[self.elvis]['S'][key][1:]

        figname = os.path.join(figspath, 'HK_%s_%s.eps' %
                               (key, self.inputs['test']))
        HKtools.doHKSinglePlot(time, HKdata[key], key, ylabel='V',
                               HKlims=HKlims, filename=figname,
                               fontsize=14)
        HKfigs.append(figname)

    self.report.add_FigsTable(HKfigs, Ncols=3, figswidth=4,
                              caption='Selected HK parameters during test.')


def save_CDP(self, cdpobj):
    cdpobj.savehardcopy()
    cdpobj.savetopickle()


def get_checkstats_T(self):
    """" """

    Xindices = copy.deepcopy(self.dd.indices)

    if 'Quad' not in Xindices.names:
        Xindices.append(core.vIndex('Quad', vals=context.Quads))

    if 'Spot' in Xindices.names:
        Xindices.pop(Xindices.names.index('Spot'))

    valini = 0.

    newcolnames = ['chk_NPIXOFF', 'chk_NPIXSAT']
    for newcolname in newcolnames:
        self.dd.initColumn(newcolname, Xindices,
                           dtype='int32', valini=valini)

    nObs = len(Xindices.get_vals('ix'))
    CCDs = Xindices.get_vals('CCD')
    Quads = Xindices.get_vals('Quad')

    if not self.drill:

        for iObs in range(nObs):
            for jCCD, CCDk in enumerate(CCDs):
                dpath = self.dd.mx['datapath'][iObs, jCCD]
                ffits = os.path.join(dpath, '%s.fits' %
                                     self.dd.mx['File_name'][iObs, jCCD])
                ccdobj = ccd.CCD(ffits)
                vstart = self.dd.mx['vstart'][iObs][jCCD]
                vend = self.dd.mx['vend'][iObs][jCCD]

                for kQ, Quad in enumerate(Quads):

                    qdata = ccdobj.get_quad(Quad, canonical=True, extension=-1)
                    NPIXOFF = len(np.where(qdata[:, vstart:vend] == 0)[0])
                    NPIXSATUR = len(np.where(qdata[:, vstart:vend] == (2**16 - 1))[0])

                    self.dd.mx['chk_NPIXOFF'][iObs, jCCD, kQ] = NPIXOFF
                    self.dd.mx['chk_NPIXSAT'][iObs, jCCD, kQ] = NPIXSATUR


def check_metrics_T(self):
    """ """

    Xindices = self.dd.indices
    CCDs = Xindices.get_vals('CCD')
    #Quads = Xindices.get_vals('Quad')

    if self.report is not None:
        self.report.add_Section(
            keyword='check_basics', Title='SATURATIONS \& LOST PIXELS', level=1)

    # FRACTION OF SATURATED PIXELS

    if self.inputs['test'] != 'PERSIST01':

        _satur_lims = self.perflimits['SATUR']
        Ncols = self.inputs['structure']['Ncols']
        satur_lims = OrderedDict()
        for CCD in CCDs:
            satur_lims[CCD] = OrderedDict()
            for icol in range(1, Ncols + 1):
                satur_lims[CCD]['col%03i' % icol] = _satur_lims[CCD]
        arrSAT = self.dd.mx['chk_NPIXSAT'][:].copy()

        NQactive = self.ccdcalc.NAXIS1 / 2. * (self.dd.mx['vend'][:] - self.dd.mx['vstart'][:])
        NQactive = np.repeat(NQactive[..., np.newaxis], 4, axis=-1)  # extend to Quads

        FracSat = arrSAT / NQactive

        _compliance_sat = self.check_stat_perCCDandCol(FracSat, satur_lims, CCDs)

        self.addComplianceMatrix2Self(_compliance_sat, 'saturation')

        if not self.IsComplianceMatrixOK(_compliance_sat):
            self.dd.flags.add('POORQUALDATA')
            self.dd.flags.add('SATURATION')
        if self.log is not None:
            self.addComplianceMatrix2Log(
                _compliance_sat, label='COMPLIANCE SATURATION FRACTION')
        if self.report is not None:
            self.addComplianceMatrix2Report(
                _compliance_sat, label='COMPLIANCE SATURATION FRACTION',
                caption='Fraction (over 1) of CCD image saturated.'
            )

    # MISSING PIXELS

    NPIX_LOST_mx = self.dd.mx['chk_NPIXOFF'][:].copy()
    NPIX_LOST = np.sum(NPIX_LOST_mx, axis=(1, 2))
    # NPIX_LOST[:] = 1 # TESTS
    ObsID = self.dd.mx['ObsID'][:].copy()

    NPIX_LOST_summ = 'LOST PIXELS: '

    if np.any(NPIX_LOST > 0):
        for iObs, Obs in enumerate(ObsID):
            _NPIX = NPIX_LOST[iObs]
            if _NPIX > 0:
                NPIX_LOST_summ += '%i (OBSID:%i), ' % (_NPIX, Obs)

        NPIX_LOST_summ = NPIX_LOST_summ[0:-2]

        self.dd.flags.add('LOSTPIXELS')

    else:

        NPIX_LOST_summ += ' None'

    if self.log is not None:
        self.log.info(NPIX_LOST_summ)

    if self.report is not None:
        self.report.add_Text(NPIX_LOST_summ)
