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
from vison.pipe import lib as pilib
import datetime
import os
import string as st
import astropy
# END IMPORT

#def generate_Explog_old(struct,defaults,elvis='6.0.0',date=pilib.dtobj_default):
#    """ """
#    
#    
#    Nscriptcols = struct['Ncols']
#    
#    columns = ELtools.columnlist[elvis]
#    
#    explog = ELtools.iniExplog(elvis)
#    
#    ixObsID = 1000
#    
#    for iscrcol in range(1,Nscriptcols+1):
#        scriptcol = struct['col%i' % iscrcol]
#        N = scriptcol['N']
#        inputkeys = [key for key in scriptcol.keys() if key != 'frames']
#        
#        rowdict = {}
#        
#        for subixrow in range(N):
#        
#            for ecol in columns:
#                rowdict[ecol] = defaults[ecol]
#        
#            for key in inputkeys:
#                rowdict[key] = scriptcol[key]
#            
#            for ixCCD in range(1,4):
#                
#                dmy = date.strftime('%d%m%y')
#                hms = date.strftime('%H%M%S')
#            
#                rowdict['ObsID'] = ixObsID
#                rowdict['File_name'] = 'EUC_%i_%sD_%sT_ROE1_CCD%i' % (ixObsID,dmy,hms,ixCCD)
#                rowdict['DATE'] = '%sD%sT' % (dmy,hms)
#                rowdict['CCD'] = 'CCD%i' % ixCCD
#                
#                explog.add_row(vals=[rowdict[key] for key in columns])
#                
#                date = date + datetime.timedelta(seconds=90)
#            
#            ixObsID += 1
#                    
#            
#    return explog


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
        assert isinstance(explog,astropy.table.table.Table)
        
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



def generate_HK(explog,defaults,datapath='',elvis='6.0.0'):
    """ """
    
    HKkeys = HKtools.allHK_keys[elvis]
    
    Nobs = len(explog['ObsID'])
    
    
    doneObsids = []
    
    for ixobs,obsid in enumerate(explog['ObsID']):
        
        if obsid in doneObsids: continue # to avoid duplications 
                                     # (each CCD has an entry in explog, so 3 entries per OBSID)
        
        idate = explog['DATE'][ixobs]
        
        idtobj = pilib.get_dtobj(idate)
        
        if ixobs < Nobs-1:
            ip1dtobj = pilib.get_dtobj(explog['DATE'][ixobs+1])
            dt = (ip1dtobj-idtobj).seconds
        else:
            dt = 90
        
        HKfilef = 'HK_%s_%s_ROE1.txt' % (obsid,idate)
        
        HKfilef = os.path.join(datapath,HKfilef)
        
        HKfile = HKtools.iniHK_QFM(elvis)
        
        for sec in range(dt):
            iidtobj = idtobj + datetime.timedelta(seconds=sec)
            
            iTimeStamp = iidtobj.strftime('%H:%M:%S')
            
            rowdict = {}
            
            for HKkey in HKkeys:
                rowdict[HKkey] = defaults[HKkey]
                
            rowdict['TimeStamp'] = iTimeStamp
            
            HKfile.add_row(vals=[rowdict[key] for key in HKkeys])
        
        
        HKfile.write(HKfilef,format='ascii',overwrite=True,delimiter='\t')
        
        doneObsids.append(obsid)
        

#def generate_FITS(wavelength,explog,datapath='',elvis='6.0.0'):
#    """ """
#    
#    NAXIS1,NAXIS2 = 4238,4132
#    
#    maxexptime = explog['Exptime'].max()
#    flatlevel = 200.
#    biaslevel = 2000.
#    pflux = 1.E4 # adu
#    pfwhm = 9./12. # pixels
#    
#    FilterID = pilib.get_FW_ID(wavelength)
#    
#    mirror_nom = polib.mirror_nom[FilterID]
#    
#    waivedkeys = ['File_name','Flsh-Rdout_e_time','C.Inj-Rdout_e_time',
#                  'Wavelength']
#    
#    for ixobs,obsid in enumerate(explog['ObsID']):
#        
#        
#        idate = explog['DATE'][ixobs]
#        iCCD = explog['CCD'][ixobs]
#        iexptime = explog['Exptime'][ixobs]
#        iMirr_pos = explog['Mirr_pos'][ixobs]
#        
#        #if iexptime == 0.: continue # TESTS
#        
#        idtobj = pilib.get_dtobj(idate)
#        
#        dmy = idtobj.strftime('%d%m%y')
#        HMS = idtobj.strftime('%H%M%S')
#
#        FITSf = 'EUC_%s_%sD_%sT_ROE1_%s.fits' % \
#            (obsid,dmy,HMS,iCCD)
#        
#        FITSf = os.path.join(datapath,FITSf)
#        
#        ccdobj = ccd.CCD()
#        
#        img = np.zeros(shape=(NAXIS1,NAXIS2),dtype='float32')
#        
#        ccdobj.add_extension(data=None)
#        ccdobj.add_extension(data=img,label='ROE1_%s' % iCCD)
#        
#        ccdobj.simadd_flatilum(levels=dict(E=flatlevel*iexptime/maxexptime*1.,
#                                           F=flatlevel*iexptime/maxexptime*1.1,
#                                           G=flatlevel*iexptime/maxexptime*1.2,
#                                           H=flatlevel*iexptime/maxexptime*1.3))
#        
#        if iexptime > 0:
#            
#            ipflux = pflux * iexptime/maxexptime
#            ipfwhm = pfwhm * (1.+((iMirr_pos-mirror_nom)/0.2)**2.)
#            
#            ccdobj.simadd_points(ipflux,ipfwhm,CCDID=iCCD,dx=0,dy=0)
#        
#        
#        ccdobj.simadd_poisson()
#        
#        ccdobj.simadd_bias(levels=dict(E=biaslevel*1.,
#                                           F=biaslevel*1.1,
#                                           G=biaslevel*1.2,
#                                           H=biaslevel*1.3))
#        ccdobj.simadd_ron()
#        
#        ccdobj.extensions[-1].data = np.round(ccdobj.extensions[-1].data).astype('int32')
#
#        ccdobj.extensions[-1].header['WAVELENG'] = explog['Wavelength'][ixobs]
#        
#        for key in ELtools.columnlist[elvis]:
#            if key not in waivedkeys:
#                ccdobj.extensions[-1].header[key] = explog[key][ixobs]
#        
#        ccdobj.writeto(FITSf,clobber=True,unsigned16bit=True)
        
    

