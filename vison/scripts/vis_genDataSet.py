"""

Development: Creating Calibration Campaign Fake Data-sets.

:History:

Created on Tue Sep 05 16:07:00 2017

:autor: 
Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import os

from vison.datamodel import EXPLOGtools as ELtools
from vison.datamodel import generator as gen
#from vison.datamodel import scriptic as sc
from vison.support import context
from vison.campaign import campaign
#from vison.point import lib as polib

import datetime

# END IMPORT


def genExpLog(toGen, explogf, equipment, elvis=context.elvis, CHAMBER=None):
    """ """

    OBSID0 = 1000

    logdefaults = {'egse_ver': elvis, 'con_file': 'vis_roe_config_cotsqm_273_vn.txt',
                   'fl_rdout': 0, 'ci_rdout': 0,
                   'fpga_ver': '2AC',
                   'R1C1_TT': -153., 'R1C2_TT': -153., 'R1C3_TT': -153.,
                   'R1C1_TB': -153., 'R1C2_TB': -153., 'R1C3_TB': -153., }

    test_sequence = campaign.generate_test_sequence(
        equipment, toGen, elvis=elvis, CHAMBER=CHAMBER)

    tests = list(test_sequence.keys())

    explog = None

    structtest0 = test_sequence[tests[0]]
    explog = gen.generate_Explog(
        structtest0,
        logdefaults,
        elvis=elvis,
        explog=explog,
        OBSID0=OBSID0,
        date=date0,
        CHAMBER=CHAMBER)

    for test in tests[1:]:
        structtest = test_sequence[test]
        explog = gen.generate_Explog(
            structtest, logdefaults, elvis=elvis, explog=explog, CHAMBER=CHAMBER)

    # WRITING EXPOSURE LOG

    explog.write(explogf, format='ascii', overwrite=True, delimiter='\t')

    return explog


def datasetGenerator(TestsSelector, doGenExplog, doGenHK, doGenFITS, outpath, elvis,
                     CHAMBER, Nrows=0):
    """ """

    equipment = dict(operator='raf',
                     sn_ccd1='CCD1TEST',
                     sn_ccd2='CCD2TEST',
                     sn_ccd3='CCD3TEST',
                     sn_roe='ROETEST',
                     sn_rpsu='RPSUTEST')

    datekey = date0.strftime('%d%m%y')
    explogf = os.path.join(outpath, 'EXP_LOG_%s.txt' % datekey)

    if doGenExplog:

        print('Generating EXPOSURE LOG...')

        explog = genExpLog(TestsSelector, explogf, equipment, elvis, CHAMBER)

    else:
        explog = ELtools.loadExpLog(explogf, elvis=elvis)

    if doGenHK:

        print('Generating HK files...')

#        HKvals = {'HK_OD_Top_CCD1':'OD_T1_V','HK_OD_Bottom_CCD1':'OD_B1_V',
#                  'HK_OD_Top_CCD2':'OD_T2_V','HK_OD_Bottom_CCD2':'OD_B2_V','HK_OD_Top_CCD3':'OD_T3_V',
#                  'HK_OD_Bottom_CCD3':'OD_B3_V','HK_IG1_Top_CCD1':'IG1_T1_V','HK_IG1_Bottom_CCD1':'IG1_B1_V',
#                  'HK_IG1_Top_CCD2':'IG1_T2_V','HK_IG1_Bottom_CCD2':'IG1_B2_V',
#                  'HK_IG1_Top_CCD3':'IG1_T3_V','HK_IG1_Bottom_CCD3':'IG1_B3_V',
#                  'HK_temp_top_CCD1':'R1CCD1TT','HK_temp_bottom_CCD1':'R1CCD1TB',
#                  'HK_temp_top_CCD2':'R1CCD2TT','HK_temp_bottom_CCD2':'R1CCD2TB',
#                  'HK_temp_top_CCD3':'R1CCD3TT','HK_temp_bottom_CCD3':'R1CCD3TB',
#                  'HK_RD_top':'RD_T_V','HK_RD_bot':'RD_B_V',
#                  'HK_IG2_top':'IG2_T_V','HK_IG2_bot':'IG2_B_V',
#                  'HK_IDH':'IDH_V','HK_IDL':'IDL_V','HK_DD_bias':24.,'HK_OG_bias':1.,'HK_1.5V_ROE':1.5,'HK_VCCD_ROE':32.,
#                  'HK_5VA_pos_ROE':5.,'HK_5V_ref_ROE':5.,'HK_10VA_ROE':10.,'HK_5.2V_neg_ROE':-5.2,'HK_3V_neg_ROE':-3.2,
#                  'HK_VRclk_ROE':10.2,'HK_VRClk_Lo_ROE':0.3,'HK_3.3V_DIG_RPSU':6.6,'HK_I3.3V_DIG_RPSU':3.3,
#                  'HK_1.5V_DIG_RPSU':3.4,'HK_I1.5V_DIG_RPSU':3.3,'HK_28V_Pri_RPSU':28.,'HK_I28V_RPSU':3.3,
#                  'HK_VAN_pos_RPSU':9.9,'HK_I+VAN_RPSU':3.3,'HK_VAN_neg_RPSU':9.9,'HK_I-VAN_RPSU':3.3,
#                  'HK_VCLK_RPSU':19.7,'HK_IVCLK_RPSU':3.3,'HK_VCCD_RPSU':50.1,'HK_IVCCD_RPSU':0.13,'HK_Temp1_RPSU':50.,
#                  'HK_Temp2_RPSU':50.,'HK_Video_TOP':50.,'HK_Video_BOT':50.,'HK_FPGA_TOP':50.,'HK_FPGA_BOT':50.,
#                  'HK_ID1':0.,'HK_ID2':0.,'HK_Viclk_ROE':0.}

        HKvals = {
            'CCD1_OD_T': 'OD_1_T',
            'CCD2_OD_T': 'OD_2_T',
            'CCD3_OD_T': 'OD_3_T',
            'COMM_RD_T': 'RD_T',
            'CCD2_IG1_T': 'IG1_2_T',
            'CCD3_IG1_T': 'IG1_3_T',
            'CCD1_TEMP_T': 'R1C1_TT',
            'CCD2_TEMP_T': 'R1C2_TT',
            'CCD3_TEMP_T': 'R1C3_TT',
            'CCD1_IG1_T': 'IG1_1_T',
            'COMM_IG2_T': 'IG2_T',
            'VID_PCB_TEMP_T': 10.,
            'FPGA_PCB_TEMP_T': 10.,
            'CCD1_OD_B': 'OD_1_B',
            'CCD2_OD_B': 'OD_2_B',
            'CCD3_OD_B': 'OD_3_B',
            'COMM_RD_B': 'RD_B',
            'CCD2_IG1_B': 'IG1_2_B',
            'CCD3_IG1_B': 'IG1_3_B',
            'CCD1_TEMP_B': 'R1C1_TB',
            'CCD2_TEMP_B': 'R1C2_TB',
            'CCD3_TEMP_B': 'R1C3_TB',
            'CCD1_IG1_B': 'IG1_1_B',
            'COMM_IG2_B': 'IG2_B',
            'VID_PCB_TEMP_B': 10.,
            'FPGA_PCB_TEMP_B': 10.,
            'FPGA_BIAS_DD': 24.,
            'FPGA_BIAS_OG': 1.,
            'FPGA_BIAS_ID1': 0.,
            'FPGA_BIAS_ID2': 0.,
            'FPGA_BIAS_ID_T': 'IDL',
            'FPGA_BIAS_ID_B': 'IDL',
            'FPGA_VRCLK_V': 10.2,
            'FPGA_10VA_P_V': 10.,
            'FPGA_VICLK_V': 0.,
            'FPGA_5VA_P_V': 5.,
            'FPGA_5VREF_V': 5.,
            'FPGA_VCCD_V': 32.,
            'FPGA_1V5VD_P_V': 1.5,
            'FPGA_3V2_N_V': -3.2,
            'FPGA_5VA_N_V': -5.2,
            'RPSU_VCCD_V': 50.1,
            'RPSU_VCLK_V': 19.7,
            'RPSU_VAN_P_V': 9.9,
            'RPSU_VAN_N_V': 9.9,
            'RPSU_3V3VD_V': 6.6,
            'RPSU_1V5VD_V': 3.4,
            'RPSU_28V_PRI_V': 28,
            'RPSU_TEMP1': 10.,
            'RPSU_VCCD_I': 0.13,
            'RPSU_VCLK_I': 19.7,
            'RPSU_VAN_P_I': 3.3,
            'RPSU_VAN_N_I': 3.3,
            'RPSU_3V3VD_I': 3.3,
            'RPSU_1V5VD_I': 3.3,
            'RPSU_28V_PRI_I': 3.3,
            'RPSU_TEMP_2': 10.,
            'stTMRErrFlg': 0,
            'hkTMRErrFlg': 0,
            'spwTmTOFlg': 0,
            'CDPUClkSt': 0,
            'fpgaSpwErr': 0,
            'V3v3ProtCnt': 0,
            'V5ProtCnt': 0,
            'VccdErrFlg': 0,
            'VclkErrFlg': 0,
            'VanErrFlg': 0,
            'stRamErr': 0,
            'hkRamErr': 0,
            'ADC_BSY_ERR_CNT': 0,
            'SPW_STATUS_REG': 12192}

        if Nrows > 0:
            sexplog = explog[0:Nrows]
        else:
            sexplog = explog
        gen.generate_HK(sexplog, HKvals, datapath=outpath, elvis=elvis)

    if doGenFITS:

        # tests
        # selix = [0,113,314,3988] # bias=0,flat=314,chinj=113,psf=3988
        #selix = [0,113,314,3988]

        print('Generating FITS files... ')

        if Nrows > 0:
            sexplog = explog[0:Nrows]
        else:
            sexplog = explog

        gen.generate_FITS_fromExpLog(sexplog, outpath, elvis, CHAMBER)


if __name__ == '__main__':

    doGenExplog = True
    doGenHK = True
    doGenFITS = True
    Nrows = None
    elvis = '6.3.0'
    CHAMBER = 'A_JUN18'

    #date0 = pilib.dtobj_default
    date0 = datetime.datetime(2018, 1, 19, 7, 0, 0)  # early riser

    # TestsSelector = dict(BIAS01=1,DARK01=1,CHINJ01=1,CHINJ02=1,
    #                  FLAT01=1,FLAT02=1,PTC01=1,PTC02WAVE=0,PTC02TEMP=0,NL01=1,
    #                  PSF01=1,PSF02=1,
    #                  TP01=1,TP02=1,
    #                  PERSIST01=1,FOCUS00=1)

    TestsSelector = dict(BIAS01=0, DARK01=0, CHINJ01=0, CHINJ02=0,
                         FLAT01=0, FLAT02=0, PTC01=0, PTC02WAVE=0, PTC02TEMP=0, NL01=0,
                         PSF01=0, PSF02=0,
                         TP01=0, TP02=0,
                         PERSIST01=0, FOCUS00=1)

    outpath = os.path.join('TEST_DATA', date0.strftime('%d_%b_%y'))

    if not os.path.exists(outpath):
        os.system('mkdir %s' % outpath)

    datasetGenerator(TestsSelector, doGenExplog, doGenHK, doGenFITS, outpath, elvis,
                     CHAMBER, Nrows=Nrows)
