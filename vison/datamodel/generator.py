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
import datetime
import os
import string as st
import astropy as astpy
from time import time
import numpy as np

# END IMPORT


Quads = ['E','F','G','H']
rdouttime = 71 # seconds, integer


def _update_fromscript(rowdict,scriptcol):
    """ """
    
    elog2sc = ELtools.script_keys_cross
    elog2sckeys = elog2sc.keys()
    
    for key in rowdict:
        if key in elog2sckeys:
            rowdict[key] = scriptcol[elog2sc[key]]
    
    phkeys = ['iphi1','iphi2','iphi3','iphi4']
    N_P_high = st.join(['I%i' % (i+1,) for i in range(4) if scriptcol[phkeys[i]]],'')    
    rowdict['N_P_high'] = N_P_high
    
    Trappump = '%i-%s' % (scriptcol['tpump'],scriptcol['tpump_mode'])
    rowdict['Trappump'] = Trappump
    
    waveix = int(scriptcol['wavelength'][-1])
    rowdict['Wavelength'] = waveix
    
    volt_cross = dict(IDL_V='IDL',IDH_V='IDH',
                          IG1_T1_V='IG1_1_T',IG1_T2_V='IG1_2_T',IG1_T3_V='IG1_3_T',
                          IG1_B1_V='IG1_1_B',IG1_B2_V='IG1_2_B',IG1_B3_V='IG1_3_B',
                          IG2_T_V='IG2_T',IG2_B_V='IG2_B',
                          OD_T1_V='OD_1_T',OD_T2_V='OD_2_T',OD_T3_V='OD_3_T',
                          OD_B1_V='OD_1_T',OD_B2_V='OD_2_T',OD_B3_V='OD_3_T',
                          RD_T_V='RD_T',RD_B_V='RD_B')
    
    for key in volt_cross.keys():
        rowdict[key] = scriptcol[volt_cross[key]]/1.E3
    
    return rowdict
    
def generate_Explog(scrdict,defaults,elvis='6.0.0',explog=None,OBSID0=1000,
                        date=pilib.dtobj_default):
    """ 
    
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
        
        datelast = HKtools.parseDTstr(explog['DATE'][-1])
        
        exptsec = explog['Exptime'][-1]/1.e3
        
        date = datelast + datetime.timedelta(seconds=80.+exptsec)
    
    
    
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
            rowdict['EGSE_ver'] = elvis
            rowdict['SPW_clk'] =  0
            rowdict['CalScrpt'] = 'Uknown'
                   
            exptsec = scriptcol['exptime']/1.e3
            
            for ixCCD in range(1,4):
                
                dmy = date.strftime('%d%m%y')
                hms = date.strftime('%H%M%S')
            
                rowdict['ObsID'] = ixObsID
                rowdict['File_name'] = 'EUC_%i_%sD_%sT_ROE1_CCD%i' % (ixObsID,dmy,hms,ixCCD)
                rowdict['DATE'] = '%sD%sT' % (dmy,hms)
                rowdict['CCD'] = 'CCD%i' % ixCCD

                       
                explog.add_row(vals=[rowdict[key] for key in expcolkeys])
                
            date = date + datetime.timedelta(seconds=80. + exptsec)
            
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
        

def generate_HK(explog,vals,datapath='',elvis='6.0.0'):
    """ """
    
        
    doneObsids = []
    
    t0 = time()
    
    for ixobs,obsid in enumerate(explog['ObsID']):
        
        if obsid in doneObsids: continue # to avoid duplications 
                                         # (each CCD has an entry in explog, so 3 entries per OBSID)
        
        idate = explog['DATE'][ixobs]
        
        idtobj = pilib.get_dtobj(idate)
        
        HKfilef = 'HK_%s_%s_ROE1.txt' % (obsid,idate)
        
        HKfilef = os.path.join(datapath,HKfilef)
        
        HKfile = HKtools.iniHK_QFM(elvis,length=int(rdouttime))

        TimeStamp = np.array([idtobj+datetime.timedelta(seconds=sec) for sec in np.arange(int(rdouttime))])        
        
        HKfile['TimeStamp'] = TimeStamp
        
        HKfile = _fill_HKcols(HKfile,explog[ixobs],vals)
        
        
        HKfile.write(HKfilef,format='ascii',overwrite=True,delimiter='\t')
        
        doneObsids.append(obsid)
        
    t1 = time()
    dtmin = (t1-t0)/60.
    nobs = len(doneObsids)
    print '%.3f minutes in generating %i HK files' % (dtmin,nobs)

    return None   


def _add_ron_window_round(ccdobj,vstart,vend):
    """ """
    
    ccdobj.simadd_ron()
    ccdobj.extensions[-1].data = np.round(ccdobj.extensions[-1].data).astype('int32')    
    
    if vstart != 0 or vend != 2066:
        ccdobj.sim_window(vstart,vend)
    
    return ccdobj

gen_bias_levels = dict(E=2300.,F=2400.,G=2500.,H=2600.)

expt_FWC_flat = dict(Filter1=2.E3,
                Filter2=2.E3,
                Filter3=2.E3,
                Filter4=2.E3,
                Filter5=2.E3,
                Filter6=200.*1.E3)

expt_FWC_point = dict(Filter1=2.E3,
                Filter2=2.E3,
                Filter3=2.E3,
                Filter4=2.E3,
                Filter5=2.E3,
                Filter6=200.*1.E3)

fwhm_nom = 9./12.

def IMG_bias_gen(ccdobj,ELdict):
    """ """
    
    vstart = ELdict['Vstart']
    vend = ELdict['Vend']

    ccdobj.simadd_bias(levels=gen_bias_levels)
    
    ccdobj = _add_ron_window_round(ccdobj,vstart,vend)

    ccdobj.extensions[-1].data = np.round(ccdobj.extensions[-1].data).astype('int32')

    return ccdobj

def IMG_flat_gen(ccdobj,ELdict):
    """ """

    vstart = ELdict['Vstart']
    vend = ELdict['Vend']
    
    wavelength = ELdict['Wavelength']
    exptime = ELdict['Exptime']
    
    tsatur = expt_FWC_flat['Filter%i' % wavelength]
    
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
    
    vstart = ELdict['Vstart']
    vend = ELdict['Vend']    
    
    chinj = ELdict['Chrg_inj']
    tpump = ELdict['Trappump']
    
    stCCD = ELdict['CCD']
    iCCD = int(stCCD[-1])
    
    noff = ELdict['Off_cycl']
    non = ELdict['On_cycle']

    IG1 = ELdict['IG1_T%i_V' % iCCD] # don't care if IG1_B is different
    IG2 = ELdict['IG2_T_V'] 
    IDL = ELdict['IDL_V']
    
    doInject = ((chinj == 1) or (tpump ==1)) and (IDL < IG1+inj_threshold)
    
    if doInject:
        
        injlevel = 2000. - max(0,IG1-IG2)/0.5*2000.
        
        injlevels = dict(E=injlevel,F=injlevel,G=injlevel,H=injlevel)
                              
        ccdobj.simadd_flatilum(levels=injlevels)
        
        ccdobj.simadd_injection(levels=injlevels,on=non,off=noff)


    ccdobj.simadd_bias(levels=gen_bias_levels) # add bias    

    ccdobj = _add_ron_window_round(ccdobj,vstart,vend)

    ccdobj.extensions[-1].data = np.round(ccdobj.extensions[-1].data).astype('int32')
    

    return ccdobj
    
def IMG_point_gen(ccdobj,ELdict):
    """ """

    vstart = ELdict['Vstart']
    vend = ELdict['Vend']
    
    wavelength = ELdict['Wavelength']
    exptime = ELdict['Exptime']
    mirror = ELdict['Mirr_pos']
    iCCD = ELdict['CCD']
    
    mirror_nom = polib.mirror_nom['Filter%i' % wavelength]    
    tsatur = expt_FWC_flat['Filter%i' % wavelength]
    
    fluence = 2.*2.**16 * exptime / tsatur
    
    fwhm = fwhm_nom * (1.+((mirror-mirror_nom)/0.2)**2.)
    
    ccdobj.simadd_points(fluence,fwhm,CCDID=iCCD,dx=0,dy=0)

    ccdobj.simadd_bias(levels=gen_bias_levels) # add bias    

    ccdobj = _add_ron_window_round(ccdobj,vstart,vend)

    ccdobj.extensions[-1].data = np.round(ccdobj.extensions[-1].data).astype('int32')

    
    return ccdobj


def generate_FITS(ELdict,funct,filename='',elvis='6.0.0'):
    """ """
    
    NAXIS1,NAXIS2 = 4238,4132
    
    
    waivedkeys = ['File_name','Flsh-Rdout_e_time','C.Inj-Rdout_e_time',
                  'Wavelength']
    CCD = ELdict['CCD']
    
    
    # CCD object initialisation
    
    ccdobj = ccd.CCD()
    img = np.zeros(shape=(NAXIS1,NAXIS2),dtype='float32')        
    ccdobj.add_extension(data=None)
    ccdobj.add_extension(data=img,label='ROE1_%s' % CCD)
    
    
    ccdobj = funct(ccdobj,ELdict)
    
    ccdobj.extensions[-1].header['WAVELENG'] = ELdict['Wavelength']
    
    for key in ELtools.columnlist[elvis]:
        if key not in waivedkeys:
            ccdobj.extensions[-1].header[key] = ELdict[key]
    
    
    if filename != '':
        ccdobj.writeto(filename,clobber=True,unsigned16bit=True)
        return None
    else:
        return ccdobj



def generate_FITS_fromExpLog(explog,datapath,elvis='6.0.0'):
    """ """    
    
    IMGgens = dict(BIAS=IMG_bias_gen,FLAT=IMG_flat_gen,
                          CHINJ=IMG_chinj_gen,POINT=IMG_point_gen)
    
    t0 = time()
    
    Nfiles = len(explog)
    
    for ixrow,obsid in enumerate(explog['ObsID']):
        
        row = explog[ixrow]

        #date = row['DATE']
        #CCD = row['CCD']
        test = row['TEST']
        
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




