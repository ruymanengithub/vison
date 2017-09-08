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

from vison.datamodel import scriptic as sc
from vison.point import FOCUS00,PSF0X
from vison.dark import BIAS01,DARK01
from vison.flat import NL01, PTC0X, FLAT0X
from vison.inject import CHINJ01,CHINJ02
from vison.pump import TP01, TP02
from vison.other import PERSIST01 as PER01
from vison.point import lib as polib

from vison.ogse.ogse import tFWC_flat,tFWC_point

import datetime

# END IMPORT

def f_write_script(struct,filename,outpath,elvis): 
    
    script = sc.Script(structure=struct,elvis=elvis)
    script.build_cargo()
    script.write(filename)    
    os.system('mv %s %s/' % (filename,outpath))   
    


def scwriter(toWrite,outpath,equipment,elvis='6.0.0'):
    """ """
    
    
    datetag = (datetime.datetime.now()).strftime('%d%b%y')
    
    
    operator = equipment['opertor']
    sn_ccd1 = equipment['sn_ccd1']
    sn_ccd2 = equipment['sn_ccd2']
    sn_ccd3 = equipment['sn_ccd3']
    sn_roe =  equipment['sn_roe']
    sn_rpsu = equipment['sn_rpsu']
        
    
    if not os.path.exists(outpath):
        os.system('mkdir %s' % outpath)
    
    # BIAS
    
    Nbias01 = 25
    fileBIAS01 = 'vis_CalCamp_BIAS01_%s_v%s.xlsx' % (datetag,elvis)
    diffBIAS01 = dict(pos_cal_mirror=polib.mirror_nom['Filter1'],sn_ccd1=sn_ccd1,
                      sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
                      sn_rpsu=sn_rpsu,operator=operator)
    
    if toWrite['BIAS01']: 
         structBIAS01 = BIAS01.build_BIAS01_scriptdict(Nbias01,
                                        diffvalues=diffBIAS01,elvis=elvis)
         f_write_script(structBIAS01,fileBIAS01,outpath,elvis)
    
    
    # DARKS
    
    Ndark01 = 4
    exptime_dark01 = 565*1.E3 # ms
    fileDARK01 = 'vis_CalCamp_DARK01_%s_v%s.xlsx' % (datetag,elvis)
    diffDARK01 = dict(pos_cal_mirror=polib.mirror_nom['Filter1'],sn_ccd1=sn_ccd1,
                      sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
                      sn_rpsu=sn_rpsu,operator=operator)
    
    if toWrite['DARK01']: 
        structDARK01 = DARK01.build_DARK01_scriptdict(Ndark01,exptime_dark01,
                                    diffvalues=diffDARK01,elvis=elvis)
        f_write_script(structDARK01,fileDARK01,outpath,elvis)


    
    # CHARGE INJECTION
    
    # CHINJ01
    
    IDL = 1100
    IDH = 1800
    IG1s = [2000,6000]
    toi_chinj01 = 500
    id_delays = [toi_chinj01*3,toi_chinj01*2]
    
    fileCHINJ01 = 'vis_CalCamp_CHINJ01_%s_v%s.xlsx' % (datetag,elvis)
    diffCHINJ01 = dict(pos_cal_mirror=polib.mirror_nom['Filter4'],sn_ccd1=sn_ccd1,
                      sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
                      sn_rpsu=sn_rpsu,operator=operator)
    
    if toWrite['CHINJ01']: 
        structCHINJ01 = CHINJ01.build_CHINJ01_scriptdict(IDL,IDH,IG1s,id_delays,
                                toi_chinj01,diffvalues=diffCHINJ01,elvis=elvis)
        f_write_script(structCHINJ01,fileCHINJ01,outpath,elvis)
    

    # CHINJ02
    
    IDLs = [13000,16000]
    IDH = 1800
    #IG1s = [2000,6000]
    toi_chinj02 = 500
    id_delays = [toi_chinj02*3,toi_chinj02*2]
    fileCHINJ02 = 'vis_CalCamp_CHINJ02_%s_v%s.xlsx' % (datetag,elvis)
    diffCHINJ02 = dict(pos_cal_mirror=polib.mirror_nom['Filter4'],sn_ccd1=sn_ccd1,
                      sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
                      sn_rpsu=sn_rpsu,operator=operator)
    
    if toWrite['CHINJ02']: 
        structCHINJ02 = CHINJ02.build_CHINJ02_scriptdict(IDLs,IDH,id_delays,toi_chinj02,
                                    diffvalues=diffCHINJ02,elvis=elvis)
        f_write_script(structCHINJ02,fileCHINJ02,outpath,elvis)     

    
    
    # FLATS
    
    t_dummy = np.array([25.,50.,75])/100.
                      
    exptimes_FLAT0X = dict(nm590=t_dummy * tFWC_flat['nm590'],
                           nm640=t_dummy * tFWC_flat['nm640'],
                           nm730=t_dummy * tFWC_flat['nm730'],
                           nm800=t_dummy * tFWC_flat['nm800'],
                           nm880=t_dummy * tFWC_flat['nm880'])
    
    # FLAT-01

    exptimesF01 = exptimes_FLAT0X['nm800'] # ms
    
    fileFLAT01 = 'vis_CalCamp_FLAT01_%s_v%s.xlsx' % (datetag,elvis)
    diffFLAT01 = dict(sn_ccd1=sn_ccd1,
                      sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
                      sn_rpsu=sn_rpsu,operator=operator,
                      test='FLAT01_800')
    
    if toWrite['FLAT01']: 

        structFLAT01 = FLAT0X.build_FLAT0X_scriptdict(exptimesF01,
                        diffvalues=diffFLAT01,elvis=elvis)
        f_write_script(structFLAT01,fileFLAT01,outpath,elvis)
    
    
    # FLAT-02
    
    wavesFLAT02 = [590,640,730,880]
    

    diffFLAT02 = dict(sn_ccd1=sn_ccd1,
                      sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
                      sn_rpsu=sn_rpsu,operator=operator)
    
    if toWrite['FLAT02']: 
    
        for iw, wave in enumerate(wavesFLAT02):
            
            iexptimes = exptimes_FLAT0X['nm%i' % wave]
            
            ifileFLAT02 = 'vis_CalCamp_FLAT02_%inm_%s_v%s.xlsx' % (wave,datetag,elvis)
            
            itestkey = 'FLAT02_%i' % wave
        
            istructFLAT02 = FLAT0X.build_FLAT0X_scriptdict(iexptimes,wavelength=wave,
                                testkey=itestkey,diffvalues=diffFLAT02,elvis=elvis)
            
            f_write_script(istructFLAT02,ifileFLAT02,outpath,elvis)
    
        
    # PTC
    
    
    
    # PTC-01
    
    filePTC01 = 'vis_CalCamp_PTC01_%s_v%s.xlsx' % (datetag,elvis)
    diffPTC01 = dict(test='PTC01',sn_ccd1=sn_ccd1,
                      sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
                      sn_rpsu=sn_rpsu,operator=operator)
        
    # 5%, 10%, 20%, 30%, 50%, 70%, 80%, 90%, 100%, 110%, 120%
    exptsPTC01 = np.array([5.,10.,20.,30.,50.,70.,80.,90.,100.,110.,120.])/100.*tFWC_flat['nm800'] # ms
    frsPTC01 = [10,10,10,10,10,10,10,10,4,4,4]
    
    if toWrite['PTC01']: 
    
        structPTC01 = PTC0X.build_PTC0X_scriptdict(exptsPTC01,frsPTC01,
                                wavelength=800,diffvalues=diffPTC01,elvis=elvis)
        f_write_script(structPTC01,filePTC01,outpath,elvis)    


    
    # PTC-02 - wavelength
    
    wavesPTC02w = [590,640,730,880]
    
    diffPTC02w = dict(sn_ccd1=sn_ccd1,
                      sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
                      sn_rpsu=sn_rpsu,operator=operator)
    
    if toWrite['PTC02WAVE']: 
    
        for iw, wave in enumerate(wavesPTC02w):
            
            # 10%, 30%, 50%, 70%, 80%, 90% x FWC. 4 frames per fluence.
            
            exptsPTC02w = np.array([10.,30.,50.,70.,80.,90.])/100.*tFWC_flat['nm%i' % wave]
            frsPTC02w = [4,4,4,4,4,4]
            
            ifilePTC02w = 'vis_CalCamp_PT0C2_%inm_%s_v%s.xlsx' % (wave,datetag,elvis)
            
            diffPTC02w['test'] = 'PTC02_%i' % wave
        
            istructPTC02w = PTC0X.build_PTC0X_scriptdict(exptsPTC02w,frsPTC02w,
                                wavelength=wave,diffvalues=diffPTC02w,elvis=elvis)
            
            f_write_script(istructPTC02w,ifilePTC02w,outpath,elvis)
    
    
   # PTC-02 - Temp.

    wavePTC02T = 800
    TempsPTC02T = [150.,156.]
    
    diffPTC02T = dict(sn_ccd1=sn_ccd1,
                      sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
                      sn_rpsu=sn_rpsu,operator=operator)

    # 10%, 30%, 50%, 70%, 80%, 90% x FWC. 4 frames per fluence.
    exptsPTC02T = np.array([10.,30.,50.,70.,80.,90.])/100.*tFWC_flat['nm%i' % wavePTC02T]
    frsPTC02T = [4,4,4,4,4,4]

    if toWrite['PTC02TEMP']: 
   
        for it,T in enumerate(TempsPTC02T):    
            
            ifilePTC02T = 'vis_CalCamp_PT0C2_T%i_%s_v%s.xlsx' % (T,datetag,elvis)
            
            diffPTC02T['test'] = 'PTC02_%iK' % T
        
            istructPTC02T = PTC0X.build_PTC0X_scriptdict(exptsPTC02T,frsPTC02T,
                                wavelength=wavePTC02T,diffvalues=diffPTC02T,elvis=elvis)
            
            f_write_script(istructPTC02T,ifilePTC02T,outpath,elvis)   
   
    
    # NL-01

    fileNL01 = 'vis_CalCamp_NL01_%s_v%s.xlsx' % (datetag,elvis)
    diffNL01 = dict(test='NL01',sn_ccd1=sn_ccd1,
                      sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
                      sn_rpsu=sn_rpsu,operator=operator)
    
    # 5 frames per fluence: 1%, 2%, 3%, 5%, 10%, 20%,30%, 50%,70%,80%,85%,90%,95%
    exptsNL01 = np.array([5.,10.,20.,30.,50.,70.,80.,90.,100.,110.,120.])/100.*tFWC_flat['ND'] # ms
    exptinterNL01 = 0.5 * tFWC_flat['ND']
    frsNL01 = np.ones(11,dtype='int32')*5

    if toWrite['NL01']: 
    
        structNL01 = NL01.build_NL01_scriptdict(exptsNL01,exptinterNL01,frsNL01,wavelength=0,
                            diffvalues=diffNL01,elvis=elvis)
        f_write_script(structNL01,fileNL01,outpath,elvis)    


    # PSF
    
    
    # PSF01-800nm
    
    wavesPSF01w = [590,640,730,800,880]
    
    diffPSF01w = dict(sn_ccd1=sn_ccd1,
                      sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
                      sn_rpsu=sn_rpsu,operator=operator)

    if toWrite['PSF01']: 

        for iw, wave in enumerate(wavesPSF01w):
            
            # 10%, 30%, 50%, 70%, 80%, 90% x FWC. 4 frames per fluence.
            
            exptsPSF01w = np.array([5.,25.,50.,75.,90.])/100.*tFWC_point['nm%i' % wave]
            frsPSF01w = [20,15,10,4,3]
            
            ifilePSF01w = 'vis_CalCamp_PSF01_%inm_%s_v%s.xlsx' % (wave,datetag,elvis)
            
            diffPSF01w['test'] = 'PSF01_%i' % wave
        
            istructPSF01w = PSF0X.build_PSF0X_scriptdict(exptsPSF01w,frsPSF01w,
                                    wavelength=wave,diffvalues=diffPSF01w,elvis=elvis)
            
            f_write_script(istructPSF01w,ifilePSF01w,outpath,elvis)     
        
        
    # PSF02
    
    wPSF02 = 800
    temps_PSF02 = [150,156]
    
    
    exptsPSF02 = np.array([5.,25.,50.,75.,90.])/100.*tFWC_point['nm%i' % wPSF02]
    frsPSF02 = [20,15,10,4,3]
    
    diffPSF02 = dict(sn_ccd1=sn_ccd1,
                      sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
                      sn_rpsu=sn_rpsu,operator=operator)
    
    
    if toWrite['PSF02']: 

    
        for it,temp in enumerate(temps_PSF02):
            
            
            ifilePSF02 = 'vis_CalCamp_PSF02_%iK_%s_v%s.xlsx' % (temp,datetag,elvis)
            
            diffPSF02['test'] = 'PSF02_%iK' % temp
        
            istructPSF02 = PSF0X.build_PSF0X_scriptdict(exptsPSF02,frsPSF02,
                        wavelength=wPSF02,diffvalues=diffPSF02,elvis=elvis)
            
            f_write_script(istructPSF02,ifilePSF02,outpath,elvis)     
    
    
    # TRAP-PUMPING

    # TP01
    
    
    TOI_TPv = [200,1000,2000,4000,8000]
    toi_chinjTP01 = 250 # quick injection
    id_delays_TP01 = np.array([3.,2.])*toi_chinjTP01
                           
    fileTP01 = 'vis_CalCamp_TP01_%s_v%s.xlsx' % (datetag,elvis)
    diffTP01 = dict(sn_ccd1=sn_ccd1,
                      sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
                      sn_rpsu=sn_rpsu,operator=operator,toi_chinj=toi_chinjTP01)

    if toWrite['TP01']: 
   
        structTP01 = TP01.build_TP01_scriptdict(Nshuffles_V=5000,TOI_TPv=TOI_TPv,
                             id_delays=id_delays_TP01,
                                diffvalues=diffTP01,elvis=elvis)
        
        f_write_script(structTP01,fileTP01,outpath,elvis)
    
    
    
    # TP02
    
    dwell_hv = [0.,4.75,14.3,28.6] # us
    toi_chinjTP02 = 250 # quick injection
    id_delays_TP02 = np.array([3.,2.])*toi_chinjTP02
          
    fileTP02 = 'vis_CalCamp_TP02_%s_v%s.xlsx' % (datetag,elvis)
    diffTP02 = dict(sn_ccd1=sn_ccd1,
                      sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
                      sn_rpsu=sn_rpsu,operator=operator,toi_chinj=toi_chinjTP02)

    if toWrite['TP02']: 
    
        structTP02 = TP02.build_TP02_scriptdict(Nshuffles_H=5000,dwell_hv=dwell_hv,
                        id_delays=id_delays_TP02,diffvalues=diffTP02,elvis=elvis)
        
        f_write_script(structTP02,fileTP02,outpath,elvis)    
    
    
    # OTHER
    
    # PERSIST
    
    exptPER01_SATUR = 50*tFWC_point['nm800']   # ms
    exptPER01_LATEN = 565.*1.E3 # ms
    
    diffPER01 = dict(sn_ccd1=sn_ccd1,
                      sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
                      sn_rpsu=sn_rpsu,operator=operator)
    
    filePER01 = 'vis_CalCamp_PERSIST01_%s_v%s.xlsx' % (datetag,elvis)
    
    if toWrite['PERSIST01']: 
    
        structPER01 = PER01.build_PER01_scriptdict(exptPER01_SATUR,exptPER01_LATEN,
                                diffvalues=diffPER01,elvis=elvis)
        
        f_write_script(structPER01,filePER01,outpath,elvis)        
    
    # FOCUS
    
    wavesFOCUS00w = [590,640,730,800,880]
    
    diffFOCUS00w = dict(sn_ccd1=sn_ccd1,
                      sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
                      sn_rpsu=sn_rpsu,operator=operator)
    

    if toWrite['FOCUS00']: 

        for iw, wave in enumerate(wavesFOCUS00w):
            
            
            iexptimeF00 = 60./100.*tFWC_point['nm%i' % wave]
            
            ifileFOCUS00w = 'vis_CalCamp_FOCUS00_%inm_%s_v%s.xlsx' % (wave,datetag,elvis)
            
            diffFOCUS00w['test'] = 'FOCUS00_%i' % wave
        
            istructFOCUS00w = FOCUS00.build_FOCUS00_scriptdict(wave,iexptimeF00,
                                    diffvalues=diffFOCUS00w,elvis=elvis)
            
            f_write_script(istructFOCUS00w,ifileFOCUS00w,outpath,elvis)         


if __name__ =='__main__':
    
    
    elvis = '6.1.0'
    
    
    outpath = 'CAL_scripts'
    
    equipment = dict(operator = 'raf',
    sn_ccd1 = 'CCD1TEST',
    sn_ccd2 = 'CCD2TEST',
    sn_ccd3 = 'CCD3TEST',
    sn_roe= 'ROETEST',
    sn_rpsu = 'RPSUTEST')
    
    toWrite = dict(BIAS01=1,DARK01=1,CHINJ01=1,CHINJ02=1,
                      FLAT01=1,FLAT02=1,PTC01=1,PTC02WAVE=1,PTC02TEMP=1,NL01=1,
                      PSF01=1,PSF02=1,
                      TP01=1,TP02=1,
                      PERSIST01=1,FOCUS00=1)
    
    scwriter(toWrite,outpath,equipment,elvis)
    
