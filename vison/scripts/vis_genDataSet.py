"""

EUCLID-VIS Ground Calibration Campaign

Development: Creating Calibration Campaign Fake Data-set

Created on Tue Sep 05 16:07:00 2017

:autor: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import os

from vison.datamodel import EXPLOGtools as ELtools
from vison.datamodel import generator as gen
from vison.datamodel import scriptic as sc
from vison.pipe import lib as pilib
from vison.point import FOCUS00,PSF0X
from vison.dark import BIAS01,DARK01
from vison.flat import NL01, PTC0X, FLAT0X
from vison.inject import CHINJ01,CHINJ02
from vison.pump import TP01, TP02
from vison.other import PERSIST01 as PER01

from vison.point import lib as polib

import datetime

# END IMPORT

    

def genExpLog(doGen,explogf,elvis):
    """ """
    
    
    OBSID0 = 1000
    operator = 'raf'
    sn_ccd1 = 'CCD1TEST'
    sn_ccd2 = 'CCD2TEST'
    sn_ccd3 = 'CCD3TEST'
    sn_roe= 'ROETEST'
    sn_rpsu = 'RPSUTEST'
    
    
    #if not os.path.exists(outpath):
    #    os.system('mkdir %s' % outpath)
    
    logdefaults = {'Lab_ver':'15.0.0','Con_file':'vis_roe_config_cotsqm_273_vn.txt',
        'CnvStart':0,'Flsh-Rdout_e_time':0,'C.Inj-Rdout_e_time':0,
        'FPGA_ver':'2AC','Chmb_pre':1.e-6,
        'R1CCD1TT':-153.,'R1CCD2TT':-153.,'R1CCD3TT':-153.,
        'R1CCD1TB':-153.,'R1CCD2TB':-153.,'R1CCD3TB':-153.,}
    
    # BIAS
    
    Nbias01 = 25
    #fileBIAS01 = 'vis_CalCamp_BIAS01_%s_v%s.xlsx' % (datetag,elvis)
    diffBIAS01 = dict(pos_cal_mirror=polib.mirror_nom['Filter1'],sn_ccd1=sn_ccd1,
                      sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
                      sn_rpsu=sn_rpsu,operator=operator)
    
    if doGen['BIAS01']:
        
        print 'BIAS01...'
        
        structBIAS01 = BIAS01.build_BIAS01_scriptdict(Nbias01,
                                    diffvalues=diffBIAS01,elvis=elvis)
        explog = gen.generate_Explog(structBIAS01,logdefaults,elvis=elvis,explog=None,OBSID0=OBSID0,
                        date=date0)
    
    # DARKS
    
    Ndark01 = 4
    exptime_dark01 = 565*1.E3 # ms
    #fileDARK01 = 'vis_CalCamp_DARK01_%s_v%s.xlsx' % (datetag,elvis)
    diffDARK01 = dict(pos_cal_mirror=polib.mirror_nom['Filter1'],sn_ccd1=sn_ccd1,
                      sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
                      sn_rpsu=sn_rpsu,operator=operator)
    
    if doGen['DARK01']: 
        
        print 'DARK01...'
        
        structDARK01 = DARK01.build_DARK01_scriptdict(Ndark01,exptime_dark01,
                    diffvalues=diffDARK01,elvis=elvis)
        
        explog = gen.generate_Explog(structDARK01,logdefaults,elvis=elvis,explog=explog)

    
    # CHARGE INJECTION
    
    # CHINJ01
    
    IDL = 1100
    IDH = 1800
    IG1s = [2000,6000]
    toi_chinj01 = 500
    id_delays = [toi_chinj01*3,toi_chinj01*2]
    
    #fileCHINJ01 = 'vis_CalCamp_CHINJ01_%s_v%s.xlsx' % (datetag,elvis)
    diffCHINJ01 = dict(pos_cal_mirror=polib.mirror_nom['Filter4'],sn_ccd1=sn_ccd1,
                      sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
                      sn_rpsu=sn_rpsu,operator=operator)
    
    if doGen['CHINJ01']: 
        
        print 'CHINJ01...'
        
        structCHINJ01 = CHINJ01.build_CHINJ01_scriptdict(IDL,IDH,IG1s,id_delays,
                    toi_chinj01,diffvalues=diffCHINJ01,elvis=elvis)
        explog = gen.generate_Explog(structCHINJ01,logdefaults,elvis=elvis,explog=explog)

    
    # CHINJ02
    
    IDLs = [13000,16000]
    IDH = 1800
    #IG1s = [2000,6000]
    toi_chinj02 = 500
    id_delays = [toi_chinj02*3,toi_chinj02*2]
    #fileCHINJ02 = 'vis_CalCamp_CHINJ02_%s_v%s.xlsx' % (datetag,elvis)
    diffCHINJ02 = dict(pos_cal_mirror=polib.mirror_nom['Filter4'],sn_ccd1=sn_ccd1,
                      sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
                      sn_rpsu=sn_rpsu,operator=operator)
    
    if doGen['CHINJ02']: 
        
        print 'CHINJ02...'
        
        structCHINJ02 = CHINJ02.build_CHINJ02_scriptdict(IDLs,IDH,id_delays,toi_chinj02,
                            diffvalues=diffCHINJ02,elvis=elvis)
        explog = gen.generate_Explog(structCHINJ02,logdefaults,elvis=elvis,explog=explog)

    
    # FLATS
    
    exptimes_FLAT0X = dict(nm590=np.array([1.,2.,3.])*1.E3,
                           nm640=np.array([1.,2.,3.])*1.E3,
                           nm730=np.array([1.,2.,3.])*1.E3,
                           nm800=np.array([1.,2.,3.])*1.E3,
                           nm880=np.array([1.,2.,3.])*1.E3)    
    
    # FLAT-01

    exptimesF01 = exptimes_FLAT0X['nm800'] # ms
    
    #fileFLAT01 = 'vis_CalCamp_FLAT01_%s_v%s.xlsx' % (datetag,elvis)
    diffFLAT01 = dict(sn_ccd1=sn_ccd1,
                      sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
                      sn_rpsu=sn_rpsu,operator=operator,
                      test='FLAT01_800')
    
    if doGen['FLAT01']: 
        
        print 'FLAT01...'

        structFLAT01 = FLAT0X.build_FLAT0X_scriptdict(exptimesF01,
                    diffvalues=diffFLAT01,elvis=elvis)
        explog = gen.generate_Explog(structFLAT01,logdefaults,elvis=elvis,explog=explog)

    
    # FLAT-02
    
    wavesFLAT02 = [590,640,730,880]
    

    diffFLAT02 = dict(sn_ccd1=sn_ccd1,
                      sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
                      sn_rpsu=sn_rpsu,operator=operator)
    
    if doGen['FLAT02']: 
    
        for iw, wave in enumerate(wavesFLAT02):
            
            
            
            iexptimes = exptimes_FLAT0X['nm%i' % wave]
            
            #ifileFLAT02 = 'vis_CalCamp_FLAT02_%inm_%s_v%s.xlsx' % (wave,datetag,elvis)
            
            itestkey = 'FLAT02_%i' % wave
        
            print '%s...' % itestkey
        
            istructFLAT02 = FLAT0X.build_FLAT0X_scriptdict(iexptimes,wavelength=wave,
                                testkey=itestkey,diffvalues=diffFLAT02,elvis=elvis)
            
            explog = gen.generate_Explog(istructFLAT02,logdefaults,elvis=elvis,explog=explog)

        
    # PTC
    
    tFWC_flat = dict(nm590=2.E3,
                nm640=2.E3,
                nm730=2.E3,
                nm800=2.E3,
                nm880=2.E3,
                ND=200.*1.E3)
    
    
    # PTC-01
    
    #filePTC01 = 'vis_CalCamp_PTC01_%s_v%s.xlsx' % (datetag,elvis)
    diffPTC01 = dict(test='PTC01',sn_ccd1=sn_ccd1,
                      sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
                      sn_rpsu=sn_rpsu,operator=operator)
        
    # 5%, 10%, 20%, 30%, 50%, 70%, 80%, 90%, 100%, 110%, 120%
    exptsPTC01 = np.array([5.,10.,20.,30.,50.,70.,80.,90.,100.,110.,120.])/100.*tFWC_flat['nm800'] # ms
    frsPTC01 = [10,10,10,10,10,10,10,10,4,4,4]
    
    if doGen['PTC01']: 
        
        print 'PTC01...'
    
        structPTC01 = PTC0X.build_PTC0X_scriptdict(exptsPTC01,frsPTC01,
                            wavelength=800,diffvalues=diffPTC01,elvis=elvis)
        explog = gen.generate_Explog(structPTC01,logdefaults,elvis=elvis,explog=explog)

    
    # PTC-02 - wavelength
    
    wavesPTC02w = [590,640,730,880]
    
    diffPTC02w = dict(sn_ccd1=sn_ccd1,
                      sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
                      sn_rpsu=sn_rpsu,operator=operator)
    
    if doGen['PTC02WAVE']: 
    
        for iw, wave in enumerate(wavesPTC02w):
            
            # 10%, 30%, 50%, 70%, 80%, 90% x FWC. 4 frames per fluence.
            
            exptsPTC02w = np.array([10.,30.,50.,70.,80.,90.])/100.*tFWC_flat['nm%i' % wave]
            frsPTC02w = [4,4,4,4,4,4]
            
            #ifilePTC02w = 'vis_CalCamp_PT0C2_%inm_%s_v%s.xlsx' % (wave,datetag,elvis)
            
            itestkey = 'PTC02_%i' % wave
            
            diffPTC02w['test'] = itestkey
            
            print '%s...' % itestkey
        
            istructPTC02w = PTC0X.build_PTC0X_scriptdict(exptsPTC02w,frsPTC02w,
                        wavelength=wave,diffvalues=diffPTC02w,elvis=elvis)
            
            explog = gen.generate_Explog(istructPTC02w,logdefaults,elvis=elvis,explog=explog)

           
    
   # PTC-02 - Temp.

    wavePTC02T = 800
    TempsPTC02T = [150.,156.]
    
    diffPTC02T = dict(sn_ccd1=sn_ccd1,
                      sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
                      sn_rpsu=sn_rpsu,operator=operator)

    # 10%, 30%, 50%, 70%, 80%, 90% x FWC. 4 frames per fluence.
    exptsPTC02T = np.array([10.,30.,50.,70.,80.,90.])/100.*tFWC_flat['nm%i' % wavePTC02T]
    frsPTC02T = [4,4,4,4,4,4]

    if doGen['PTC02TEMP']: 
   
        for it,T in enumerate(TempsPTC02T):    
            
            #ifilePTC02T = 'vis_CalCamp_PT0C2_T%i_%s_v%s.xlsx' % (T,datetag,elvis)
            
            itestkey = 'PTC02_%iK' % T
                        
            diffPTC02T['test'] = itestkey
                      
            print '%s...' % itestkey
        
            istructPTC02T = PTC0X.build_PTC0X_scriptdict(exptsPTC02T,frsPTC02T,
                                wavelength=wavePTC02T,diffvalues=diffPTC02T,elvis=elvis)
            
            explog = gen.generate_Explog(istructPTC02T,logdefaults,elvis=elvis,explog=explog)

    
    # NL-01

    #fileNL01 = 'vis_CalCamp_NL01_%s_v%s.xlsx' % (datetag,elvis)
    diffNL01 = dict(test='NL01',sn_ccd1=sn_ccd1,
                      sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
                      sn_rpsu=sn_rpsu,operator=operator)
    
    # 5 frames per fluence: 1%, 2%, 3%, 5%, 10%, 20%,30%, 50%,70%,80%,85%,90%,95%
    exptsNL01 = np.array([5.,10.,20.,30.,50.,70.,80.,90.,100.,110.,120.])/100.*tFWC_flat['ND'] # ms
    exptinterNL01 = 0.5 * tFWC_flat['ND']
    frsNL01 = np.ones(11,dtype='int32')*5

    if doGen['NL01']: 
        
        print 'NL01...'
    
        structNL01 = NL01.build_NL01_scriptdict(exptsNL01,exptinterNL01,frsNL01,wavelength=0,
                            diffvalues=diffNL01,elvis=elvis)
        explog = gen.generate_Explog(structNL01,logdefaults,elvis=elvis,explog=explog)

    
    # PSF

    tFWC_point = dict(nm590=2.E3,
                nm640=2.E3,
                nm730=2.E3,
                nm800=2.E3,
                nm880=2.E3,
                ND=200.*1.E3)
    
    
    # PSF01-800nm
    
    wavesPSF01w = [590,640,730,800,880]
    
    diffPSF01w = dict(sn_ccd1=sn_ccd1,
                      sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
                      sn_rpsu=sn_rpsu,operator=operator)

    if doGen['PSF01']: 

        for iw, wave in enumerate(wavesPSF01w):
            
            # 10%, 30%, 50%, 70%, 80%, 90% x FWC. 4 frames per fluence.
            
            exptsPSF01w = np.array([5.,25.,50.,75.,90.])/100.*tFWC_point['nm%i' % wave]
            frsPSF01w = [20,15,10,4,3]
            
            #ifilePSF01w = 'vis_CalCamp_PSF01_%inm_%s_v%s.xlsx' % (wave,datetag,elvis)
            
            itestkey = 'PSF01_%i' % wave
            diffPSF01w['test'] = itestkey
                      
            print '%s...' % itestkey
        
            istructPSF01w = PSF0X.build_PSF0X_scriptdict(exptsPSF01w,frsPSF01w,
                                wavelength=wave,diffvalues=diffPSF01w,elvis=elvis)
            
            explog = gen.generate_Explog(istructPSF01w,logdefaults,elvis=elvis,explog=explog)

        
    # PSF02
    
    wPSF02 = 800
    temps_PSF02 = [150,156]
    
    
    exptsPSF02 = np.array([5.,25.,50.,75.,90.])/100.*tFWC_point['nm%i' % wPSF02]
    frsPSF02 = [20,15,10,4,3]
    
    diffPSF02 = dict(sn_ccd1=sn_ccd1,
                      sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
                      sn_rpsu=sn_rpsu,operator=operator)
    
    
    if doGen['PSF02']: 

    
        for it,temp in enumerate(temps_PSF02):
            
            
            #ifilePSF02 = 'vis_CalCamp_PSF02_%iK_%s_v%s.xlsx' % (temp,datetag,elvis)
            
            itestkey = 'PSF02_%iK' % temp
            
            diffPSF02['test'] = itestkey
                     
            print '%s...' % itestkey
        
            istructPSF02 = PSF0X.build_PSF0X_scriptdict(exptsPSF02,frsPSF02,
                                wavelength=wPSF02,diffvalues=diffPSF02,elvis=elvis)
            
            explog = gen.generate_Explog(istructPSF02,logdefaults,elvis=elvis,explog=explog)

            
    
    # TRAP-PUMPING

    # TP01
    
    
    TOI_TPv = [200,1000,2000,4000,8000]
    toi_chinjTP01 = 250 # quick injection
    id_delays_TP01 = np.array([3.,2.])*toi_chinjTP01
                           
    #fileTP01 = 'vis_CalCamp_TP01_%s_v%s.xlsx' % (datetag,elvis)
    diffTP01 = dict(sn_ccd1=sn_ccd1,
                      sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
                      sn_rpsu=sn_rpsu,operator=operator,toi_chinj=toi_chinjTP01)

    if doGen['TP01']: 
        
        print 'TP01...'
   
        structTP01 = TP01.build_TP01_scriptdict(Nshuffles_V=5000,TOI_TPv=TOI_TPv,id_delays=id_delays_TP01,
                                   diffvalues=diffTP01,elvis=elvis)
        
        explog = gen.generate_Explog(structTP01,logdefaults,elvis=elvis,explog=explog)

    
    # TP02
    
    dwell_hv = [0.,4.75,14.3,28.6] # us
    toi_chinjTP02 = 250 # quick injection
    id_delays_TP02 = np.array([3.,2.])*toi_chinjTP02
          
    #fileTP02 = 'vis_CalCamp_TP02_%s_v%s.xlsx' % (datetag,elvis)
    diffTP02 = dict(sn_ccd1=sn_ccd1,
                      sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
                      sn_rpsu=sn_rpsu,operator=operator,toi_chinj=toi_chinjTP02)

    if doGen['TP02']: 
        
        print 'TP02...'
    
        structTP02 = TP02.build_TP02_scriptdict(Nshuffles_H=5000,dwell_hv=dwell_hv,
            id_delays=id_delays_TP02,diffvalues=diffTP02,elvis=elvis)
        
        explog = gen.generate_Explog(structTP02,logdefaults,elvis=elvis,explog=explog)

    
    # OTHER
    
    # PERSIST
    
    exptPER01_SATUR = 15*1.E3   # ms
    exptPER01_LATEN = 565.*1.E3 # ms
    
    diffPER01 = dict(sn_ccd1=sn_ccd1,
                      sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
                      sn_rpsu=sn_rpsu,operator=operator)
    
    #filePER01 = 'vis_CalCamp_PERSIST01_%s_v%s.xlsx' % (datetag,elvis)
    
    if doGen['PERSIST01']: 
        
        print 'PERSIST01...'
    
        structPER01 = PER01.build_PER01_scriptdict(exptPER01_SATUR,exptPER01_LATEN,
                    diffvalues=diffPER01,elvis=elvis)
        
        explog = gen.generate_Explog(structPER01,logdefaults,elvis=elvis,explog=explog)

   
    # FOCUS
    
    wavesFOCUS00w = [590,640,730,800,880]
    
    diffFOCUS00w = dict(sn_ccd1=sn_ccd1,
                      sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
                      sn_rpsu=sn_rpsu,operator=operator)
    

    if doGen['FOCUS00']: 

        for iw, wave in enumerate(wavesFOCUS00w):
            
            
            iexptimeF00 = 60./100.*tFWC_point['nm%i' % wave]
            
            
            itestkey = 'FOCUS00_%i' % wave
            
            diffFOCUS00w['test'] = itestkey
                        
            print '%s...' % itestkey
        
            istructFOCUS00w = FOCUS00.build_FOCUS00_scriptdict(wave,iexptimeF00,
                            diffvalues=diffFOCUS00w,elvis=elvis)
            
            explog = gen.generate_Explog(istructFOCUS00w,logdefaults,elvis=elvis,explog=explog)

     
    # WRITING EXPOSURE LOG
    
    explog.write(explogf,format='ascii',overwrite=True,delimiter='\t')
    
    return explog

    
    
def datasetGenerator(TestsSelector,doGenExplog,doGenHK,doGenFITS,outpath,elvis,
                     Nrows=0):
    """ """
    
    
    datekey = date0.strftime('%d%m%y')    
    explogf = os.path.join(outpath,'EXP_LOG_%s.txt' % datekey)
    
    
    if doGenExplog:
        explog = genExpLog(TestsSelector,explogf,elvis)
    else:
        explog = ELtools.loadExpLog(explogf,elvis=elvis)
    
    
    
    if doGenHK:
        
        HKvals = {'HK_OD_Top_CCD1':'OD_T1_V','HK_OD_Bottom_CCD1':'OD_B1_V',
                  'HK_OD_Top_CCD2':'OD_T2_V','HK_OD_Bottom_CCD2':'OD_B2_V','HK_OD_Top_CCD3':'OD_T3_V',
                  'HK_OD_Bottom_CCD3':'OD_B3_V','HK_IG1_Top_CCD1':'IG1_T1_V','HK_IG1_Bottom_CCD1':'IG1_B1_V',
                  'HK_IG1_Top_CCD2':'IG1_T2_V','HK_IG1_Bottom_CCD2':'IG1_B2_V',
                  'HK_IG1_Top_CCD3':'IG1_T3_V','HK_IG1_Bottom_CCD3':'IG1_B3_V',
                  'HK_temp_top_CCD1':'R1CCD1TT','HK_temp_bottom_CCD1':'R1CCD1TB',
                  'HK_temp_top_CCD2':'R1CCD2TT','HK_temp_bottom_CCD2':'R1CCD2TB',
                  'HK_temp_top_CCD3':'R1CCD3TT','HK_temp_bottom_CCD3':'R1CCD3TB',
                  'HK_RD_top':'RD_T_V','HK_RD_bot':'RD_B_V',
                  'HK_IG2_top':'IG2_T_V','HK_IG2_bot':'IG2_B_V',
                  'HK_IDH':'IDH_V','HK_IDL':'IDL_V','HK_DD_bias':24.,'HK_OG_bias':1.,'HK_1.5V_ROE':1.5,'HK_VCCD_ROE':32.,
                  'HK_5VA_pos_ROE':5.,'HK_5V_ref_ROE':5.,'HK_10VA_ROE':10.,'HK_5.2V_neg_ROE':-5.2,'HK_3V_neg_ROE':-3.2,
                  'HK_VRclk_ROE':10.2,'HK_VRClk_Lo_ROE':0.3,'HK_3.3V_DIG_RPSU':6.6,'HK_I3.3V_DIG_RPSU':3.3,
                  'HK_1.5V_DIG_RPSU':3.4,'HK_I1.5V_DIG_RPSU':3.3,'HK_28V_Pri_RPSU':28.,'HK_I28V_RPSU':3.3,
                  'HK_VAN_pos_RPSU':9.9,'HK_I+VAN_RPSU':3.3,'HK_VAN_neg_RPSU':9.9,'HK_I-VAN_RPSU':3.3,
                  'HK_VCLK_RPSU':19.7,'HK_IVCLK_RPSU':3.3,'HK_VCCD_RPSU':50.1,'HK_IVCCD_RPSU':0.13,'HK_Temp1_RPSU':50.,
                  'HK_Temp2_RPSU':50.,'HK_Video_TOP':50.,'HK_Video_BOT':50.,'HK_FPGA_TOP':50.,'HK_FPGA_BOT':50.,
                  'HK_ID1':0.,'HK_ID2':0.,'HK_Viclk_ROE':0.}
        
        if Nrows > 0:
            sexplog = explog[0:Nrows]
        else:
            sexplog = explog
        gen.generate_HK(sexplog,HKvals,datapath=outpath,elvis=elvis)
    
    
    if doGenFITS:
        
        #tests
        #selix = [0,113,314,3988] # bias=0,flat=314,chinj=113,psf=3988        
        #selix = [0,113,314,3988]        
        
        if Nrows > 0:
            sexplog = explog[0:Nrows]
        else:
            sexplog = explog
        
        gen.generate_FITS_fromExpLog(sexplog,outpath,elvis)
     

if __name__ =='__main__':
    
    
    
    doGenExplog = True
    doGenHK = True
    doGenFITS = True
    Nrows = 0
    
    date0 = pilib.dtobj_default 
    elvis = '6.1.0'
    
    TestsSelector = dict(BIAS01=1,DARK01=0,CHINJ01=1,CHINJ02=0,
                      FLAT01=1,FLAT02=0,PTC01=0,PTC02WAVE=0,PTC02TEMP=0,NL01=1,
                      PSF01=1,PSF02=0,
                      TP01=1,TP02=0,
                      PERSIST01=0,FOCUS00=1)
    
    outpath = os.path.join('TEST_DATA',date0.strftime('%d_%b_%y'))
    
    datasetGenerator(TestsSelector,doGenExplog,doGenHK,doGenFITS,outpath,elvis,
                     Nrows=Nrows)
