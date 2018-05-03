# -*- coding: utf-8 -*-
"""

Support Script with Variables and Functions used for FM-Calib. analysis

Created on Wed Nov 30 11:11:27 2016


:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk
"""

# IMPORT STUFF
from pdb import set_trace as stop

import copy
from glob import glob
import os
import numpy as np
import string as st
import datetime

from vison import __version__
from vison.datamodel import core
from vison.datamodel import EXPLOGtools as ELtools
from vison.datamodel import HKtools
from vison.support import context
#from vison.support import files
#from vison.support import vistime
# END IMPORT


#elvis = '6.5.X'
#Quads =['E','F','G','H']
# RON = 4.5 # e-
# gain = 3.5 # e-/ADU

#ccdobj = CCD()
#prescan = ccdobj.prescan
#overscan = ccdobj.overscan
#imgheight = ccdobj.NAXIS2/2
#quad_width = ccdobj.NAXIS1/2
#imgwidth = quad_width - prescan - overscan

# sumwell = dict(fwd_bas=[9.425,4.825],fwd_e2v=[7.475,6.55],
#               rwd_bas_v=[6.5,0.175],rwd_bas_s=[6.5,0.175],
#               rwd_bas_vs=[6.5,0.175])

def extract_date_explogf(explogf):
    """ """
    items = st.split(explogf, '_')
    stdate = items[2]
    idate = datetime.datetime.strptime(stdate, '%d%m%y')
    return idate


def sortbydateexplogfs(explogfs):
    """Sorts exposure logs according to date"""

    assert isinstance(explogfs, list)

    dates = []
    for explogf in explogfs:
        baseexplogf = os.path.basename(explogf)
        baseexplogf = os.path.splitext(baseexplogf)[0]
        idate = extract_date_explogf(baseexplogf)
        dates.append(idate)
    return (np.array(explogfs)[np.argsort(dates)]).tolist()


def loadexplogs(explogfs, elvis=context.elvis, addpedigree=False, datapath=None):
    """loads in memory an explog (text-file) or list of explogs (text-files)."""

    if isinstance(explogfs, str):
        explog = ELtools.loadExpLog(explogfs, elvis=elvis)

    elif isinstance(explogfs, list):
        expLogList = []
        for explogf in explogfs:
            expLogList.append(ELtools.loadExpLog(explogf, elvis=elvis))
        explog = ELtools.mergeExpLogs(expLogList, addpedigree)

    # add datapath(s)
    # The datapath becomes another column in DataDict. This helps dealing
    # with tests that run overnight and for which the input data is in several
    # date-folders.

    if datapath is not None:

        if isinstance(datapath, list):
            assert len(datapath) == len(explogfs)

            longestdatapathname = max([len(item) for item in datapath])
            explog['datapath'] = np.zeros(
                len(explog), dtype='S%i' % longestdatapathname)
            explognumber = explog['explognumber']
            for idata in range(len(datapath)):
                explog['datapath'][explognumber == idata] = datapath[idata]
        else:
            explog['datapath'] = np.array([datapath]*len(explog))

    return explog


def check_test_structure(explog, structure, CCDs=[1, 2, 3], selbool=True, wavedkeys=[]):
    """Checks whether a selected number of exposures is consistent with the 
    expected acquisition test structure.


    """

    from vison.datamodel.generator import _update_fromscript

    isconsistent = True

    Ncols = structure['Ncols']

    elogkeys = explog.keys()

    dummyrowdict = dict([(key, None) for key in elogkeys])

    if not np.any(selbool) or len(explog) == 0:

        isconsistent = False
        failedkeys = ['all']
        failedcols = np.arange(1, Ncols+1).tolist()

        report = dict(checksout=isconsistent,
                      failedkeys=failedkeys, failedcols=failedcols)

        return report

    failedkeys = []
    failedcols = []

    for iCCD in CCDs:

        cselbool = selbool & (explog['CCD'] == 'CCD%i' % iCCD)

        ix0 = 0
        for iC in range(1, Ncols+1):

            colname = 'col%i' % iC

            expectation = structure[colname]

            frames = expectation['frames']

            if frames > 0:
                ix1 = ix0+frames
            else:
                ix1 = None

            if (ix1 is not None) and (ix1-ix0 != frames):
                failedcols.append(iC)
                continue
            try:
                ixsubsel = np.where(cselbool)[0][ix0:ix1]
            except IndexError:
                failedcols.append(iC)
                continue

            rowdict = _update_fromscript(dummyrowdict, expectation)

            for key in rowdict:
                val = rowdict[key]
                if val is None or key in wavedkeys:
                    continue
                try:
                    checksout = np.all(np.isclose(explog[key][ixsubsel], val))
                except TypeError:
                    checksout = np.all(explog[key][ixsubsel] == val)
                isconsistent &= checksout
                if ~checksout:
                    failedkeys.append(key)

            ix0 += frames

    failedkeys = np.unique(failedkeys).tolist()
    failedcols = np.unique(failedcols).tolist()

    report = dict(checksout=isconsistent,
                  failedkeys=failedkeys, failedcols=failedcols)

    return report


# def oldDataDict_builder(explog,inputs,structure):
#    """ """
#
#    # Build DataDict - And Labelling
#
#    DataDict = dict(meta = dict(inputs=inputs,structure=structure))
#
#    for CCDindex in [1,2,3]:
#
#        CCDkey = 'CCD%i' % CCDindex
#
#        DataDict[CCDkey] = dict()
#
#        CCDselbool = explog['CCD'] == CCDkey
#
#        if len(np.where(CCDselbool)[0]) == 0:
#            continue
#
#        for key in explog.colnames:
#            DataDict[CCDkey][key] = explog[key][CCDselbool].data.copy()
#
#    return DataDict

def DataDict_builder(explog, inputs, structure):
    """ """

    # Build DataDict

    dd = core.DataDict()
    # Load Metadata
    dd.meta = dict(inputs=inputs, structure=structure, vison=__version__)
    # Load Exposure Log
    dd.loadExpLog(explog)

    return dd


def addHK(dd, HKKeys, elvis=context.elvis):
    """Adds HK information to a DataDict object."""

    if len(HKKeys) == 0:
        return dd

    ObsIDs = dd.mx['ObsID']().copy()
    datapaths = dd.mx['datapath'][:, 0].copy()

    HKlist = []

    for iOBS, ObsID in enumerate(ObsIDs):

        tmp = os.path.join(datapaths[iOBS], 'HK_%s_*_ROE1.txt' % ObsID)

        HKs = glob(tmp)

        if len(HKs) == 1:
            HKlist.append(HKs[0])
        elif len(HKs) > 1:
            print 'More than one HK file for ObsID %i' % ObsID
            print HKs
        elif len(HKs) == 0:
            print 'HK file for ObsID %i not found' % ObsID

    obsids, dtobjs, tdeltasec, readHKKeys, HKdata = HKtools.parseHKfiles(
        HKlist, elvis=elvis)

    HKix = copy.deepcopy(dd.mx['ObsID'].indices)

    for HKKey in HKKeys:

        pre_HKKey = 'HK_%s' % HKKey

        dd.initColumn(pre_HKKey, HKix, dtype='float32', valini=np.nan)

        ixkey = readHKKeys.index(HKKey)
        dd.mx[pre_HKKey][:] = HKdata[:, 0, ixkey]
    
    return dd


def coarsefindTestinExpLog(explog, testkey, Nframes):
    """Checks whether the 'Nframes' of test with key 'testkey' are in 'explog'.
    :param explog: astropy table, exposure log object.
    :param testkey: char, 
    :param Nframes: number of frames expected from test.
    :return: bool, test was acquired or not.

    """
    #if testkey=='BIAS01': stop()
    indices = np.where(explog['test'] == testkey)
    if len(indices[0]) == 0:
        return False

    ccds = explog['CCD'][indices]
    uccd = np.unique(ccds)
    nccds = len(uccd)
    Nobsids = len(indices[0])/nccds

    wasacquired = Nobsids == Nframes
    return wasacquired


def broadcast_todo_flags_func(inputdict,tododict):
    
    for taskname in inputdict['tasks']:
        for key in inputdict[taskname]['todo_flags']:
            inputdict[taskname]['todo_flags'][key] = False
    
    for taskname in inputdict['tasks']:
        inputdict[taskname]['todo_flags'].update(tododict)
    
    return inputdict

def broadcast_todo_flags(inputdict, docheck=False, dotest=False, doreport=False):

    assert np.array([docheck, dotest]).sum(
    ) <= 1, 'At most 1 kwd should be True!'
    if np.array([docheck, dotest]).sum() == 0:
        return inputdict
    
    _todocheck = dict(init=False,check=False,report=False)
    _todocheck['report'] = doreport
    
    if docheck:
        _todocheck.update(dict(init=True, check=True, report=True))
    elif dotest:
        _todocheck.update(dict(init=True))
        
    inputdict = broadcast_todo_flags_func(inputdict,_todocheck)

    return inputdict


