#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Created on Wed Apr 25 10:26:58 2018

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
from pdb import set_trace as stop
# END IMPORT

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


def get_test_duration(teststruct):
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