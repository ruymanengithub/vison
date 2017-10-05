#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: FLAT0X

Flat-fields acquisition / analysis script

Created on Tue Aug 29 17:32:52 2017

:author: Ruyman Azzollini
:contact: r.azzollini__at__ucl.ac.uk

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
from vison.pipe import lib as pilib
from vison.ogse import ogse
from vison.point import lib as polib
from vison.datamodel import  scriptic as sc
#from vison.pipe import FlatFielding as FFing
#from vison.support.report import Report
#from vison.support import files
from vison.datamodel import ccd
from vison.datamodel import EXPLOGtools as ELtools
from vison.datamodel import HKtools
from vison.datamodel import ccd
from vison.datamodel import generator
import datetime
from copy import deepcopy

# END IMPORT

isthere = os.path.exists


FLAT0X_commvalues = dict(program='CALCAMP',
  IPHI1=1,IPHI2=1,IPHI3=1,IPHI4=0,
  rdmode='Normal',
  flushes=7,exptime=0.,shuttr=1,
  siflush=0,
  wave=4,
  comments='')
  


def build_FLAT0X_scriptdict(exptimes,wavelength=800,testkey='FLAT0X',
                            diffvalues=dict(),elvis='6.3.0'):
    """Builds FLAT0X script structure dictionary.
    
    :param exptimes: list of ints, exposure times.
    :param wavelength: int, wavelength.
    :param testkey: char, test identifier.
    :param diffvalues: dict, opt, differential values.
    
    """
    
    FW_ID = ogse.get_FW_ID(wavelength)
    FW_IDX = int(FW_ID[-1])
    
    FLAT0X_commvalues['wavelength'] = 'Filter %i' % FW_IDX
    FLAT0X_commvalues['test'] = testkey
    
    
    FLAT0X_sdict = dict(col1=dict(frames=80,exptime=exptimes[0],comments='25pc'),
                        col2=dict(frames=60,exptime=exptimes[1],comments='50pc'),
                        col3=dict(frames=30,exptime=exptimes[2],comments='75pc'))

    Ncols = len(FLAT0X_sdict.keys())    
    FLAT0X_sdict['Ncols'] = Ncols
                
    commvalues = deepcopy(sc.script_dictionary[elvis]['defaults'])
    commvalues.update(FLAT0X_commvalues)
    
    FLAT0X_sdict = sc.update_structdict(FLAT0X_sdict,commvalues,diffvalues)
    
    return FLAT0X_sdict



def check_data(DataDict,report,inputs,log=None):
    """ 
    METACODE
    
    Checks quality of ingested data.
    
    # check common HK values are within safe / nominal margins
    # check voltages in HK match commanded voltages, within margins
    
    # f.e.ObsID, f.e.CCD, f.e.Q.:
        # measure offsets/means in pre-, img-, over-
        # measure std in pre-, img-, over-
    # assess std in pre- is within allocated margins
    # assess offsets in pre- and over- are equal, within allocated  margins
    # assess fluences are within allocated margins
    # flag saturations if there are.
    
    # plot fluence vs. time for each exptime
    # plot std-pre vs. time
    
    # issue any warnings to log
    # issue update to report

    
    """
    
    
    return DataDict, report
    

def do_indiv_flats(DataDict,report,inputs,log=None):
    """
    METACODE
    
    Preparation of data for further analysis and 
    produce flat-field for each OBSID.

    # f.e. ObsID, f.e.CCD, f.e.Q:
        # subtract offset
        # opt: [sub bias frame]
        # model 2D fluence distro in image area
        # produce average profile along rows
        # produce average profile along cols
        # save 2D model and profiles in a pick file for each OBSID-CCD
        #
        # divide by 2D model to produce indiv-flat
        # save indiv-Flat to FITS, update add filename

    # plot average profiles f. each CCD and Q (color coded by time)
    
    
    """
    
    
    return DataDict,report
    
def do_master_flat(DataDict,report,inputs,log=None):
    """ 
    METACODE
    
    Produces Master Flat-Field

    # f.e.CCD, f.e.Q:
        # stack individual flat-fields by chosen estimator
        # save Master FF to FITS
        # measure PRNU and 
    
    # report PRNU figures
    
    """
    
    return DataDict,report


def do_prdef_mask(DataDict,report,inputs,log=None):
    """
    METACODE
    
    Produces mask of defects in Photo-Response
    
    # f.e.CCD, f.e.Q:
        # produce mask of PR defects
        # save mask of PR defects
        # count dead pixels / columns 
    
    # report PR-defects stats
    
    """    



def feeder(inputs,elvis='6.3.0'):
    """ """
    
    subtasks = [('check',check_data),('indivflats',do_indiv_flats),
                ('masterflat',do_master_flat),
                ('prmask',do_prdef_mask)]
    
    
    exptimes = inputs['exptimes']
    wavelength = inputs['wavelength']
    testkey = inputs['testkey']
    if 'elvis' in inputs:
        elvis = inputs['elvis']
    if 'diffvalues' in inputs:
        diffvalues = inputs['diffvalues']
    else:
        diffvalues = {}
    
    scriptdict = build_FLAT0X_scriptdict(exptimes,wavelength,testkey,
                            diffvalues=diffvalues,elvis=elvis)

    
    inputs['structure'] = scriptdict
    inputs['subtasks'] = subtasks
    
    return inputs