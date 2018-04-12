"""

VIS Ground Calibration Campaign


Automatically Generating Calibration Campaign Scripts.

Created on Fri Sep 08 12:03:00 2017

:autor: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import os
from collections import OrderedDict
import copy

import vison
from vison.datamodel import scriptic as sc
from vison.pipe import campaign, minicampaign
from vison.support import context

import datetime
import string as st

# END IMPORT


def f_write_script(struct, filename, outpath, elvis):

    script = sc.Script(structure=struct, elvis=elvis)
    script.build_cargo()
    script.write(filename)

    stdout = os.popen('md5sum %s' % filename).read()
    xsum = st.split(stdout)[0]

    os.system('mv %s %s/' % (filename, outpath))

    return xsum


def get_dt_ro(**kwargs):
    """ """
    vstart = kwargs['vstart']
    vend = kwargs['vend']
    toi = kwargs['toi_ro']*1.e-6
    dtserial = 14.3E-6 * (51.+2048.+29.)
    nlines = vend-vstart
    dt = nlines * (4.*toi+dtserial)
    return dt


def get_dt_flush(**kwargs):
    """ """
    nlines = 2066
    toi_fl = kwargs['toi_fl'] * 1.E-6
    dt = nlines * (4.*toi_fl)
    return dt


def get_dt_chinj(**kwargs):
    """ """

    dtserial = 14.3E-6 * (51.+2048.+29.)
    vstart = 0
    vend = 2066
    nlines = vend-vstart
    toi_ch = kwargs['toi_ch']*1.e-6
    chin_dly = kwargs['chin_dly']
    dtserial *= chin_dly

    dt = nlines * (4.*toi_ch+dtserial)
    return dt


def get_dt_stpump(**kwargs):
    """ """
    s_tp_cnt = kwargs['s_tp_cnt']
    dwell_s = kwargs['dwell_s']*1.e-6
    vstart = kwargs['vstart']
    vend = kwargs['vend']
    nlines = vend-vstart

    dt = nlines * ((2.*14.3E-6+dwell_s)*s_tp_cnt)
    return dt


def get_dt_vtpump(**kwargs):
    """ """
    toi_tp = kwargs['toi_tp']*1.e-6
    v_tp_cnt = kwargs['v_tp_cnt']
    dwell_v = kwargs['dwell_v']*1.e-6
    dt = ((4.*toi_tp+dwell_v)*v_tp_cnt)

    return dt


def test_timer(teststruct):
    """Provides estimate of test duration."""

    testduration = 0.
    Ncols = teststruct['Ncols']

    for ixcol in range(1, Ncols+1):
        col = teststruct['col%i' % ixcol]
        iduration = col['exptime']

        if col['siflsh'] == 1:
            iduration += col['siflsh_p']*1.e-3

        iduration += get_dt_flush(**col)

        iduration += get_dt_ro(**col)

        if col['chinj'] == 1 or col['s_tpump'] == 1 or col['v_tpump'] == 1:
            iduration += get_dt_chinj(**col)

        if col['s_tpump'] == 1:
            iduration += get_dt_stpump(**col)
        if col['v_tpump'] == 1:
            iduration += get_dt_vtpump(**col)

        iduration *= col['frames']
        #print iduration
        testduration += iduration

    testduration /= 60.  # minutes

    return testduration


def scwriter(toWrite, test_generator, outpath, equipment, elvis=context.elvis):
    """ """

    datetag = (datetime.datetime.now()).strftime('%d%b%y')
    versiontag = vison.__version__

    checksumf = 'CHECK_SUMS_%s.txt' % datetag
    inventoryf = 'TESTS_INVENTORY_%s.txt' % datetag
    checksums = []

    if not os.path.exists(outpath):
        os.system('mkdir %s' % outpath)

    test_sequence = test_generator(equipment, toWrite, elvis=elvis)

    Nframes = 0
    duration = 0

    f1 = open(os.path.join(outpath, inventoryf), 'w')
    print >> f1, 'Scripts written on %s' % datetag
    print >> f1, 'checksumf: %s' % checksumf
    print >> f1, 'vison version: %s\n' % versiontag

    print 'WRITING SCRIPTS...'

    for test in test_sequence.keys():
        print '%s...' % test
        testobj = copy.deepcopy(test_sequence[test])
        structtest = testobj.build_scriptdict(elvis=elvis)
        testduration = test_timer(structtest)

        duration += testduration

        iNcols = structtest['Ncols']
        frameslist = [structtest['col%i' % i]['frames']
                      for i in range(1, iNcols+1)]
        iNframes = np.sum(frameslist)

        summary = '%s: %i [%.2f min] cols: %s' % (test, iNcols, testduration,
                                                  frameslist.__repr__())
        print >> f1, summary

        Nframes += iNframes

        ffile = 'vis_CalCamp_%s_%s_v%s.xlsx' % (test, datetag, elvis)
        xsum = f_write_script(structtest, ffile, outpath, elvis)
        checksums.append((ffile, xsum))

    print >> f1, '\n %i Frames Total' % Nframes
    print >> f1, '\n %.2f Minutes Total' % duration

    f1.close()

    # WRITING CHECKSUMS

    f2 = open(os.path.join(outpath, checksumf), 'w')
    for item in checksums:
        print >> f2, '%-60s\t%s' % item
    f2.close()


if __name__ == '__main__':

    elvis = '6.5.X'

    camptype = 'Full'  # Mini/Full

    outpath = 'MiniCal_data02/Full_Camp_scripts_26JAN18_E6.5.0'

    equipment = dict(operator='cpf',
                     sn_ccd1='D01',
                     sn_ccd2='D02',
                     sn_ccd3='D03',
                     sn_roe='SIM',
                     sn_rpsu='PSUs')

    toWrite_def = OrderedDict(BIAS01=0, DARK01=0, CHINJ00=0, CHINJ01=0, CHINJ02=0,
                              FLAT01=0, FLAT02=0, PTC01=0, PTC02WAVE=0, PTC02TEMP=0, NL01=0,
                              PSF01=0, PSF02=0,
                              TP00=0, TP01=0, TP02=0,
                              PERSIST01=0, FOCUS00=0)

    if camptype == 'Full':

        #toWrite = toWrite_def.copy()
        # toWrite.update(OrderedDict(BIAS01=1,DARK01=1,CHINJ01=1,CHINJ02=1,
        #              FLAT01=1,FLAT02=1,PTC01=1,PTC02WAVE=1,PTC02TEMP=1,NL01=1,
        #              PSF01=1,PSF02=1,
        #              TP01=1,TP02=1,
        #              PERSIST01=1,FOCUS00=1))

        toWrite = OrderedDict(BIAS01=1, DARK01=1, CHINJ01=1, CHINJ02=1,
                              FLAT01=1, FLAT02=1, PTC01=1, PTC02WAVE=1, PTC02TEMP=1, NL01=1,
                              PSF01=1, PSF02=1,
                              TP01=1, TP02=1,
                              PERSIST01=1, FOCUS00=1)

        test_generator = campaign.generate_test_sequence

    elif camptype == 'Mini':

        toWrite = toWrite_def.copy()
        toWrite.update(OrderedDict(BIAS01=1, DARK01=1, CHINJ00=1, TP00=1,
                                   FLAT01=1, FLAT02=1, FOCUS00=1, PSF01=1))

        test_generator = minicampaign.generate_reduced_test_sequence

    scwriter(toWrite, test_generator, outpath, equipment, elvis)
