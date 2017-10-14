#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Script to generate simulated data for pipeline testing purposes.

Created on Tue Aug 29 11:08:56 2017

:author: Ruyman Azzollini
:contact: r.azzollini__at__ucl.ac.uk

"""


# IMPORT STUFF
from pdb import set_trace as stop
from vison.datamodel import EXPLOGtools as ELtools
from vison.datamodel import HKtools
from vison.datamodel import ccd
from vison.pipe import lib as pilib
from vison.point import lib as polib
from vison.support import vistime
import datetime
import os
import string as st
import astropy as astpy
from time import time
import numpy as np

from vison.ogse import ogse

# END IMPORT

gen_bias_levels = dict(E=2300.,F=2400.,G=2500.,H=2600.)

#expt_FWC_flat = dict(Filter1=2.E3,
#                Filter2=2.E3,
#                Filter3=2.E3,
#                Filter4=2.E3,
#                Filter5=2.E3,
#                Filter6=200.*1.E3)
#
#expt_FWC_point = dict(Filter1=2.E3,
#                Filter2=2.E3,
#                Filter3=2.E3,
#                Filter4=2.E3,
#                Filter5=2.E3,
#                Filter6=200.*1.E3)
#
#fwhm_nom = 9./12. # pixels



Quads = ['E','F','G','H']
rdouttime = 71 # seconds, integer


def _update_fromscript(rowdict,scriptcol):
    """ """
    
    #elog2sc = ELtools.script_keys_cross
    #elog2sckeys = elog2sc.keys()
    
    scriptkeys = scriptcol.keys()
    
    for key in rowdict:
        if key in scriptkeys:
            rowdict[key] = scriptcol[key]
    
    phkeys = ['IPHI1','IPHI2','IPHI3','IPHI4']
    
    IPHI = st.join(['I%i' % (i+1,) for i in range(4) if scriptcol[phkeys[i]]],'')    
    rowdict['IPHI'] = IPHI
    
    #Trappump = '%i-%s' % (scriptcol['tpump'],scriptcol['tpump_mode'])
    #rowdict['Trappump'] = Trappump
    
    #waveix = int(scriptcol['wavelength'][-1])
    #rowdict['Wavelength'] = waveix
    
    #volt_cross = dict(IDL_V='IDL',IDH_V='IDH',
    #                      IG1_T1_V='IG1_1_T',IG1_T2_V='IG1_2_T',IG1_T3_V='IG1_3_T',
    #                      IG1_B1_V='IG1_1_B',IG1_B2_V='IG1_2_B',IG1_B3_V='IG1_3_B',
    #                      IG2_T_V='IG2_T',IG2_B_V='IG2_B',
    #                      OD_T1_V='OD_1_T',OD_T2_V='OD_2_T',OD_T3_V='OD_3_T',
    #                      OD_B1_V='OD_1_T',OD_B2_V='OD_2_T',OD_B3_V='OD_3_T',
    #                      RD_T_V='RD_T',RD_B_V='RD_B')
    
    #for key in volt_cross.keys():
    #    rowdict[key] = scriptcol[volt_cross[key]]/1.E3
    
    return rowdict
    
def generate_Explog(scrdict,defaults,elvis='6.3.0',explog=None,OBSID0=1000,
                        date=vistime.dtobj_default):
    """ 
    
    Generates a fake ExposureLog from a test structure dictionary.
    
    DEVELOPMENT NOTES:
        
    To be generated: (EASY)
        *ObsID, *File_name, *CCD, *ROE=ROE1, *DATE, *BUNIT=ADU,
        SPW_clk=0?,EGSE_ver=elvis,
    
    Temporal:
        SerRdDel
    
    To be provided in defaults: (EASY)   
        Lab_ver,Con_file,CnvStart,
        Flsh-Rdout_e_time,C.Inj-Rdout_e_time,
        FPGA_ver,Chmb_pre,R1CCD[1,2,3]T[T,B]
        
    To be read/parsed/processed from struct: (DIFFICULT)
        
        SerRDel?,SumWell?,
        IniSweep?,+etc.
    
    
    """
    
    
    Nscriptcols = scrdict['Ncols']
    
    expcolkeys = ELtools.columnlist[elvis]
    if explog is None:
        explog = ELtools.iniExplog(elvis)
        ixObsID = OBSID0
        
    else:
        assert isinstance(explog,astpy.table.table.Table)
        
        ixObsID = explog['ObsID'][-1]+1
        
        datelast = HKtools.parseDTstr(explog['date'][-1])
        
        #exptsec = explog['exptime'][-1]
        
        date = datelast #+ datetime.timedelta(seconds=4.+rdouttime+exptsec)
    
    
    
    for iscrcol in range(1,Nscriptcols+1):
        
        scriptcol = scrdict['col%i' % iscrcol]
        N = scriptcol['frames']
        inputkeys = [key for key in scriptcol.keys() if key != 'frames']
        
        rowdict = {}
        for eckey in expcolkeys: rowdict[eckey] = None
        
        for subixrow in range(N):
            
            rowdict.update(defaults)
            #for ecol in columns:
            #    rowdict[ecol] = defaults[ecol]
            
            rowdict = _update_fromscript(rowdict,scriptcol)
            
            for key in inputkeys:
                rowdict[key] = scriptcol[key]
            
            
            rowdict['ROE'] = 'ROE1'
            rowdict['BUNIT'] = 'ADU' 
            rowdict['egse_ver'] = elvis
            rowdict['spw_clk'] =  0
            rowdict['calscrpt'] = 'Unknown'
                   
            exptsec = scriptcol['exptime']
            
            optime = 0.
            if (scriptcol['chinj'] == 1) or (scriptcol['v_tpump'] == 1) or (scriptcol['s_tpump'] == 1):
                optime = rdouttime
            
            date = date + datetime.timedelta(seconds= 4.+rdouttime + optime + exptsec)
            
            for ixCCD in range(1,4):
                
                dmy = date.strftime('%d%m%y')
                hms = date.strftime('%H%M%S')
            
                rowdict['ObsID'] = ixObsID
                rowdict['File_name'] = 'EUC_%i_%sD_%sT_ROE1_CCD%i' % (ixObsID,dmy,hms,ixCCD)
                rowdict['date'] = '%sD%sT' % (dmy,hms)
                rowdict['CCD'] = 'CCD%i' % ixCCD

                       
                try: explog.add_row(vals=[rowdict[key] for key in expcolkeys])
                except: stop()
                
            
            ixObsID += 1
                    
            
    return explog

#def _fill_HKrow(row,vals): # obsolete
#    """ """
#    rowdict = {}
#    
#    for key in vals:
#        val = vals[key]
#        if isinstance(val,str):
#            rowdict[key] = row[val]
#        else:
#            rowdict[key] = val
#    return rowdict


def _fill_HKcols(HKfile,row,vals):
    """ """
    
    nrows = len(HKfile)
    
    for key in vals:
        val = vals[key]
        dtype = HKfile[key].dtype
        unity = np.ones(nrows,dtype=dtype)
        
        if isinstance(val,str):
            HKfile[key] = row[val] * unity
        else:
            HKfile[key] = val * unity
                   
    return HKfile

def merge_HKfiles(HKfilefs,masterHKf):
    """ """
    
    hdr = open(HKfilefs[0],'r').readlines()[0][:-1]
    
    body = []
    
    for HKfilef in HKfilefs:
        lines = open(HKfilef,'r').readlines()[1:]
        body += lines
    
    
    f = open(masterHKf,'w')
    print >> f, hdr
    for line in body:
        print >> f, line[:-1]
    f.close()
    


def generate_HK(explog,vals,datapath='',elvis='6.3.0'):
    """ """
    
    date0 = explog['date'][0]
    dtobj0 = vistime.get_dtobj(date0)
    
    masterHKf = 'HK_%s_ROE1.txt' % (dtobj0.strftime('%d-%m-%y'),)
    masterHKf = os.path.join(datapath,masterHKf)
    
    skip = False #True for some TESTS
        
    doneObsids = []
    HKfilefs = []
    
    t0 = time()
    
    
    for ixobs,obsid in enumerate(explog['ObsID']):
        
        if obsid in doneObsids: continue # to avoid duplications 
                                         # (each CCD has an entry in explog, so 3 entries per OBSID)
        
        idate = explog['date'][ixobs]
        
        idtobj = vistime.get_dtobj(idate)
        
        HKfilef = 'HK_%s_%s_ROE1.txt' % (obsid,idate)
        
        HKfilef = os.path.join(datapath,HKfilef)
        HKfilefs.append(HKfilef)
        
        if not skip:
        
            HKfile = HKtools.iniHK_QFM(elvis,length=int(rdouttime))
    
            TimeStampRaw = np.array([idtobj+datetime.timedelta(seconds=sec) for sec in np.arange(int(rdouttime))])
            TimeStamp = [item.strftime('%d-%m-%y_%H:%M:%S') for item in TimeStampRaw]
            
            HKfile['TimeStamp'] = TimeStamp
            
            HKfile = _fill_HKcols(HKfile,explog[ixobs],vals)
            
            HKfile.write(HKfilef,format='ascii',overwrite=True,delimiter='\t')
        
        doneObsids.append(obsid)

    
    #HKfilefs = HKfilefs[0:10]
    merge_HKfiles(HKfilefs,masterHKf)
    
    t1 = time()
    dtmin = (t1-t0)/60.
    nobs = len(doneObsids)
    print '%.3f minutes in generating %i HK files' % (dtmin,nobs)

    return None   


def _add_ron_window_round(ccdobj,vstart,vend):
    """ """
    
    ccdobj.simadd_ron()
    ccdobj.extensions[-1].data = np.round(ccdobj.extensions[-1].data).astype('int32')    
    
    if vstart != 1 or vend != ccd.NAXIS2/2:
        ccdobj.sim_window(vstart-1,vend)
    
    return ccdobj



def IMG_bias_gen(ccdobj,ELdict):
    """ """
    
    vstart = ELdict['vstart']
    vend = ELdict['vend']

    ccdobj.simadd_bias(levels=gen_bias_levels)
    
    ccdobj = _add_ron_window_round(ccdobj,vstart,vend)

    ccdobj.extensions[-1].data = np.round(ccdobj.extensions[-1].data).astype('int32')

    return ccdobj

def IMG_flat_gen(ccdobj,ELdict):
    """ """

    vstart = ELdict['vstart']
    vend = ELdict['vend']
    
    waveID = ELdict['wave']
    exptime = ELdict['exptime']
    
    tsatur = ogse.tFWC_flat['nm%i' % ogse.FW['F%i' % waveID]]
    
    fluence = 2.**16 * exptime / tsatur
    
    ilumlevels = dict(E=fluence,F=fluence,G=fluence,H=fluence)
    
    ccdobj.simadd_flatilum(levels=ilumlevels)
        
    ccdobj.simadd_poisson()

    ccdobj.simadd_bias(levels=gen_bias_levels) # add bias    

    ccdobj = _add_ron_window_round(ccdobj,vstart,vend)

    ccdobj.extensions[-1].data = np.round(ccdobj.extensions[-1].data).astype('int32')

    return ccdobj
    
def IMG_chinj_gen(ccdobj,ELdict,):
    """ """
    
    inj_threshold = 7.5
    
    vstart = ELdict['vstart']
    vend = ELdict['vend']    
    
    chinj = ELdict['chinj']
    
    stCCD = ELdict['CCD']
    iCCD = int(stCCD[-1])
    
    noff = ELdict['chinj_of']
    non = ELdict['chinj_on']

    IG1 = ELdict['IG1_%i_T' % iCCD] # don't care if IG1_B is different
    IG2 = ELdict['IG2_T'] 
    IDL = ELdict['IDL']
    
    doInject = (chinj == 1) and (IDL < IG1+inj_threshold)
    
    if doInject:
        
        
        injlevel = 2000. + max(0,IG2-IG1)/0.5*2000.
        
        injlevels = dict(E=injlevel,F=injlevel,G=injlevel,H=injlevel)
                              
        ccdobj.simadd_flatilum(levels=injlevels)
        
        ccdobj.simadd_injection(levels=injlevels,on=non,off=noff)


    ccdobj.simadd_bias(levels=gen_bias_levels) # add bias    

    ccdobj = _add_ron_window_round(ccdobj,vstart,vend)

    ccdobj.extensions[-1].data = np.round(ccdobj.extensions[-1].data).astype('int32')
    

    return ccdobj
    
def IMG_point_gen(ccdobj,ELdict):
    """ """

    vstart = ELdict['vstart']
    vend = ELdict['vend']
    
    waveID = ELdict['wave']
    wavenm = ogse.FW['F%i' % waveID]
    exptime = ELdict['exptime']
    mirror = ELdict['mirr_pos']
    iCCD = ELdict['CCD']
    
    mirror_nom = polib.mirror_nom['F%i' % waveID]    
    tsatur = ogse.tFWC_point['nm%i' % wavenm]
    
     
    fluence = 2.*2.**16 * exptime / tsatur
    
    fwhm = ogse.fwhm_lambda['nm%i' % wavenm] * (1.+((mirror-mirror_nom)/0.2)**2.)
    
    ccdobj.simadd_points(fluence,fwhm,CCDID=iCCD,dx=0,dy=0)

    ccdobj.simadd_bias(levels=gen_bias_levels) # add bias    

    ccdobj = _add_ron_window_round(ccdobj,vstart,vend)

    ccdobj.extensions[-1].data = np.round(ccdobj.extensions[-1].data).astype('int32')

    
    return ccdobj


def generate_FITS(ELdict,funct,filename='',elvis='6.3.0'):
    """ """
    
    NAXIS1,NAXIS2 = ccd.NAXIS1,ccd.NAXIS2
        
    waivedkeys = ['File_name','fl_rdout','ci_rdout',
                  'wave']
    
    CCD = ELdict['CCD']
    
    
    # CCD object initialisation
    
    ccdobj = ccd.CCD()
    img = np.zeros(shape=(NAXIS1,NAXIS2),dtype='float32')        
    ccdobj.add_extension(data=None)
    ccdobj.add_extension(data=img,label='ROE1_%s' % CCD)
    
    
    ccdobj = funct(ccdobj,ELdict)
    
    ccdobj.extensions[-1].header['WAVELENG'] = ELdict['wave']
    
    for key in ELtools.columnlist[elvis]:
        if key not in waivedkeys:
            ccdobj.extensions[-1].header[key] = ELdict[key]
    
    
    if filename != '':
        ccdobj.writeto(filename,clobber=True,unsigned16bit=True)
        return None
    else:
        return ccdobj



def generate_FITS_fromExpLog(explog,datapath,elvis='6.3.0'):
    """ """    
    
    IMGgens = dict(BIAS=IMG_bias_gen,FLAT=IMG_flat_gen,
                          CHINJ=IMG_chinj_gen,POINT=IMG_point_gen)
    
    t0 = time()
    
    Nfiles = len(explog)
    
    #explog = explog[0:100] # tests
    #explog = explog[3195:] # tests
    
    for ixrow,obsid in enumerate(explog['ObsID']):
        
        row = explog[ixrow]

        #date = row['DATE']
        #CCD = row['CCD']
        test = row['test']
        
        tmpFITSname = row['File_name']
        FITSname = os.path.join(datapath,'%s.fits' % tmpFITSname)
        
        isbias = np.any([key in test for key in ['BIAS','DARK']])
        if isbias: FITSgen = IMGgens['BIAS']
        
        isinj = np.any([key in test for key in ['CHINJ','TP']])
        if isinj: FITSgen = IMGgens['CHINJ']
        
        isflat = np.any([key in test for key in ['FLAT','PTC','NL']])
        if isflat: FITSgen = IMGgens['FLAT']
        
        ispoint = np.any([key in test for key in ['FOCUS','PSF','PERSIST']])
        if ispoint: FITSgen = IMGgens['POINT']
        
        generate_FITS(row,funct=FITSgen,filename=FITSname,elvis=elvis)
        
        
    t1 = time()
    dtmin = (t1-t0)/60.
    
    print '%.3f minutes in generating %i FITS files' % (dtmin,Nfiles)

    return None       




