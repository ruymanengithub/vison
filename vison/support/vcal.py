#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Module with tools for loading and applying ROE CCD voltage 
calibrations.

Created on Thu May 16 13:31:53 2019

@author: raf
"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import os
import pandas as pd
from collections import OrderedDict
import sys
import copy

# END IMPORT

DAC_keys = ['DAC Slope', 'DAC Intercept']
ADC_keys = ['ADC Slope',' ADC Intercept']

CCDs = [1,2,3]
Quads = ['E','F','G','H']


def get_board_key(Q):
    if Q in ['E','F']: return 'T'
    elif Q in ['G','H']: return 'B'


def get_Vkey_OD(CCD,Q):
    """ """
    Vkey = 'CCD%i_OD_%s_%s' % (CCD,get_board_key(Q),Q)
    return Vkey
    

def get_Vkey_RD(CCD,Q):
    """ """
    Vkey = 'CCD%i_COMM_RD_%s_%s' % (CCD,get_board_key(Q),Q)
    return Vkey
    
def get_Vkey_IG1(CCD,Q):
    """ """
    Vkey = 'CCD%i_IG1_%s' % (CCD,get_board_key(Q))
    return Vkey
    
def get_Vkey_IG2(CCD,Q):
    """ """
    Vkey = 'CCD%i_COMM_IG2_%s' % (CCD,get_board_key(Q))
    return Vkey
    
def get_Vkey_ID(CCD,Q):
    Vkey = 'CCD%i_COMM_ID_H_%s' % (CCD,get_board_key(Q))
    return Vkey


V_translators = dict(OD=get_Vkey_OD,
     RD=get_Vkey_RD,
     IG1=get_Vkey_IG1,
     IG2=get_Vkey_IG2,
     IDL=get_Vkey_ID,
     IDH=get_Vkey_ID)


CCDBIAS_ELVIS_dict = OrderedDict()
CCDBIAS_ELVIS_dict['7.5.X'] = OrderedDict(
        ccd1_od_t = [0.00723,0.296],
        ccd2_od_t = [0.007196,0.243],
        ccd3_od_t = [0.007187,0.283],
        comm_rd_t = [0.0048,0.26],
        ccd1_ig1_t = [0.0035,-0.0055],
        ccd2_ig1_t = [0.0035,-0.0084],
        ccd3_ig1_t = [0.0035,-0.0191],
        comm_ig2_t = [0.0035,0.2796],
        ccd1_od_b = [0.007199,0.281],
        ccd2_od_b = [0.00721,0.251],
        ccd3_od_b = [0.00722,0.275],
        comm_rd_b = [0.0048,0.26],
        ccd1_ig1_b = [0.0035,-0.0012],
        ccd2_ig1_b = [0.0035,-0.0013],
        ccd3_ig1_b = [0.0035,-0.0005],
        comm_ig2_b = [0.0035,0.2504],
        comm_id_h = [0.00592,0.],
        comm_id_l = [0.00592,0.]
        )

def f_V2counts_ESCRIPT(V,EVkey,version='7.5.X'):
    """ """
    p=CCDBIAS_ELVIS_dict[version][EVkey]
    return (V+p[1])/p[0]

HK_ELVIS_dict = OrderedDict()
HK_ELVIS_dict['7.5.X'] = OrderedDict(
        ccd3_od_t = [0.010904,0.],
        ccd2_od_t = [0.010904,0.],
        ccd1_od_t = [0.010904,0.],
        comm_rd_t = [0.007238,0.],
        ccd2_ig1_t = [0.005534,0.],
        ccd1_ig1_t = [0.005534,0.],
        ccd3_ig1_t = [0.005544,0.],
        comm_ig2_t = [0.005123,0.],
        ccd1_od_b = [0.010904,0.],
        ccd2_od_b = [0.010904,0.],
        ccd3_od_b = [0.010904,0.],
        comm_rd_b = [0.00726,0.],
        ccd2_ig1_b = [0.005544,0.],
        ccd3_ig1_b = [0.005507,0.],
        ccd1_ig1_b = [0.005507,0.],
        comm_ig2_b = [0.005123,0.],
        fpga_bias_id1 = [0.00591,0.],
        fpga_bias_id2 = [0.005895,0.],
        )

def f_V2counts_EHK(EV,HKkey,version='7.5.X'):
    """ """
    p = HK_ELVIS_dict[version][HKkey]
    return (EV-p[1])/p[0]



def load_VCal(calfile,ROE):
    """ """
    
    dfVCAL = pd.read_excel(calfile,sheetname='TRANSFER_FM%i' % ROE)
    
    if len(dfVCAL) == 0:
        return None
    
    VCALdict = OrderedDict()
    
    VCALdict['DAC'] = OrderedDict()
    VCALdict['ADC'] = OrderedDict()
    
    
    def _get_SlIn(df,Vkey,convk):
        slope = float(df['%s Slope' % convk].loc[Vkey])        
        icept = float(df['%s Intercept' % convk].loc[Vkey])
        return [slope, icept]
    
    
    for CCD in CCDs:
        
        CCDk = 'CCD%i' % CCD
        VCALdict['DAC'][CCDk] = OrderedDict()
        VCALdict['ADC'][CCDk] = OrderedDict()
        
        for Q in Quads:
            
            VCALdict['DAC'][CCDk][Q] = OrderedDict()
            VCALdict['ADC'][CCDk][Q] = OrderedDict()
            
            for convk in ['DAC','ADC']:
            
                Vkey_OD = get_Vkey_OD(CCD,Q)
                VCALdict[convk][CCDk][Q]['OD'] = _get_SlIn(dfVCAL,Vkey_OD,convk)
                
                Vkey_RD = get_Vkey_RD(CCD,Q)
                VCALdict[convk][CCDk][Q]['RD'] = _get_SlIn(dfVCAL,Vkey_RD,convk)
                
                Vkey_IG1 = get_Vkey_IG1(CCD,Q)
                VCALdict[convk][CCDk][Q]['IG1'] = _get_SlIn(dfVCAL,Vkey_IG1,convk)
                            
                Vkey_IG2 = get_Vkey_IG2(CCD,Q)
                VCALdict[convk][CCDk][Q]['IG2'] = _get_SlIn(dfVCAL,Vkey_IG2,convk)
                
                Vkey_ID = get_Vkey_ID(CCD,Q)
                VCALdict[convk][CCDk][Q]['ID'] = _get_SlIn(dfVCAL,Vkey_ID,convk)
    
    return VCALdict



def f_counts2V_CAL(Ecounts, VCALdict, keytuple, mode='DAC'):
    
    coeffs = VCALdict[mode][keytuple[0]][keytuple[1]][keytuple[2]]
    return Ecounts * coeffs[0] + coeffs[1]

def get_Ekey_OD(CCD,Q,mode):
    """ """
    if mode == 'SCRIPT':
        Ekey = 'OD_%s_%s' % (CCD,get_board_key(Q))
    elif mode in ['FPGA','HK']:
        Ekey = 'ccd%i_od_%s' % (CCD,get_board_key(Q).lower())
    return Ekey
    

def get_Ekey_RD(CCD,Q,mode):
    """ """
    if mode == 'SCRIPT':
        Ekey = 'RD_%s' % (get_board_key(Q),)
    elif mode in ['FPGA','HK']:
        Ekey = 'comm_rd_%s' % (get_board_key(Q).lower(),)
        
    return Ekey
    
def get_Ekey_IG1(CCD,Q,mode):
    """ """
    if mode == 'SCRIPT':
        Ekey = 'IG1_%i_%s' % (CCD,get_board_key(Q))
    elif mode in ['FPGA','HK']:
        Ekey = 'ccd%i_ig1_%s' % (CCD,get_board_key(Q).lower())
    return Ekey
    
def get_Ekey_IG2(CCD,Q,mode):
    """ """
    if mode == 'SCRIPT':
        Ekey = 'IG2_%s' % (get_board_key(Q),)
    elif mode in ['FPGA','HK']:
        Ekey = 'comm_ig2_%s' % (get_board_key(Q).lower(),)
    return Ekey
    
def get_Ekey_IDL(CCD,Q,mode):
    if mode == 'SCRIPT':
        Ekey = 'IDL'
    elif mode == 'FPGA':
        Ekey = 'comm_id_l'
    elif mode == 'HK':
        Ekey = 'fpga_bias_id2'
    return Ekey

def get_Ekey_IDH(CCD,Q,mode):
    if mode == 'SCRIPT':
        Ekey = 'IDH'
    elif mode == 'FPGA':
        Ekey = 'comm_id_h'
    elif mode == 'HK':
        Ekey = 'fpga_bias_id1'
    return Ekey

ELVIS_translators = dict(OD=get_Ekey_OD,
        RD=get_Ekey_RD,
        IG1=get_Ekey_IG1,
        IG2=get_Ekey_IG2,
        IDL=get_Ekey_IDL,
        IDH=get_Ekey_IDH)

class RoeVCal(object):
    """ """
    
    CCDs = CCDs
    Quads = Quads
    
    _ELVIS_translators = copy.deepcopy(ELVIS_translators)
    
    cal_keys = ['OD','RD','IG1','IG2','IDL','IDH']
    
    def __init__(self,calfile=None):
        """ """
        self.calfile = calfile
        self.VCaldict = dict()
        self.ROE = None
        self.Eversion = '7.5.X'
        
    def load_VCal(self,ROE):
        """ """
        self.ROE = ROE
        VCaldict = load_VCal(self.calfile, self.ROE)
        if VCaldict is None:
            raise IOError
        self.VCaldict = VCaldict
    
    
    
    def fcal_ELVIS_script(self, inV, calkey, CCD, Q):
        """Converts ELVIS input (via script ) voltages to calibrated voltages."""
        
        CCDk = 'CCD%i' % CCD
        
        ftrans_ELVIS = self._ELVIS_translators[calkey]
        
        Fkey = ftrans_ELVIS(CCD,Q,mode='FPGA')
        
        Ecounts = self.f_V2counts_ESCRIPT(inV, Fkey)
        
        keytuple = self._get_keytuple(calkey, CCDk, Q)
        
        outV = self.f_counts2V_CAL(Ecounts, keytuple, mode='DAC')
        
        return outV
    
    def get_HKkey(self,calkey, CCD, Q):
        """ """
        
        ftrans_ELVIS = self._ELVIS_translators[calkey]  
        HKkey = ftrans_ELVIS(CCD,Q,mode='HK')
        return HKkey
    
    def fcal_HK(self, inV, calkey, CCD, Q):
        """Converts ELVIS HK voltages to calibrated voltages."""
        
        CCDk = 'CCD%i' % CCD
        
        ftrans_ELVIS = self._ELVIS_translators[calkey]  
        
        HKkey = ftrans_ELVIS(CCD,Q,mode='HK')
        
        countsHK = self.f_V2counts_EHK(inV, HKkey)
        
        keytuple = self._get_keytuple(calkey, CCDk, Q)
        
        trueVHK = self.f_counts2V_CAL(countsHK, keytuple, mode='ADC')
        
        return trueVHK
    
    def _get_keytuple(self, calkey, CCDk, Q):
        """ """
        
        if calkey in ['IDL','IDH']:
            _mkey = 'ID'
        else:
            _mkey = calkey
        
        keytuple = (CCDk,Q,_mkey)
        
        return keytuple
        
    
    def f_counts2V_CAL(self,countsHK, keytuple, mode='DAC'):
        """ """
        
        return f_counts2V_CAL(countsHK, self.VCaldict, keytuple, mode)
    
    def f_V2counts_ESCRIPT(self, inV, Fkey):
        """ """
        return f_V2counts_ESCRIPT(inV, Fkey, version = self.Eversion)
    
    def f_V2counts_EHK(self, inV, HKkey):
        """ """
        return f_V2counts_EHK(inV, HKkey, version=self.Eversion)
    
    
def proto(calfile,ROE):
    """ """
    
    
    VCALdict = load_VCal(calfile,ROE)
    
    if VCALdict is None:
        print('vcal not available roe ROE %i' % ROE)
        sys.exit()
    
    script = OrderedDict(
            IDL=13,IDH=18.,
            IG1_1_T=6.,IG1_2_T=6.,IG1_3_T=6.,
            IG1_1_B=6.,IG1_2_B=6.,IG1_3_B=6.,
            IG2_T=2.,IG2_B=2.,
            OD_1_T=26.,OD_1_B=26.,
            OD_2_T=26.,OD_2_B=26.,
            OD_3_B=26.,OD_3_T=26.,
            RD_T=16.,RD_B=16.,
            comm_rd_t=16.         
            )
    
    mkeys = ['OD','RD','IG1','IG2','IDL','IDH']
    
    ELVIS_translators = dict(OD=get_Ekey_OD,
                         RD=get_Ekey_RD,
                         IG1=get_Ekey_IG1,
                         IG2=get_Ekey_IG2,
                         IDL=get_Ekey_IDL,
                         IDH=get_Ekey_IDH)
    
    print('\nCommanded vs. applied Voltages')

    for CCD in CCDs:
        CCDk = 'CCD%i' % CCD
        for Q in Quads:
            
            for mkey in mkeys:
                
                ftrans_ELVIS = ELVIS_translators[mkey]            
                Ekey = ftrans_ELVIS(CCD,Q,mode='SCRIPT')
                Fkey = ftrans_ELVIS(CCD,Q,mode='FPGA')
                
                inV = script[Ekey]
                
                Ecounts = f_V2counts_ESCRIPT(inV, Fkey, version='7.5.X')
                
                
                ftrans_CAL = V_translators[mkey]
                Vkey = ftrans_CAL(CCD,Q)
                
                if mkey in ['IDL','IDH']:
                    _mkey = 'ID'
                else:
                    _mkey = mkey
                
                keytuple = (CCDk,Q,_mkey)
                
                outV = f_counts2V_CAL(Ecounts, VCALdict, keytuple, mode='DAC')
                
                print('%s=%.1f     %s=%.2f' % (Ekey, inV, Vkey, outV))
    
    
    HK = dict(ccd3_od_t=26.115, ccd2_od_t=26.061, ccd1_od_t=26.006,
              comm_rd_t=16.09,
              ccd2_ig1_t=6.231, ccd1_ig1_t=6.22, ccd3_ig1_t=6.204,
              comm_ig2_t=2.198,
              ccd1_od_b=26.061, ccd2_od_b=26.028, ccd3_od_b=26.039,
              comm_rd_b=16.139, ccd2_ig1_b=6.253, ccd3_ig1_b=6.19,
              ccd1_ig1_b=6.206, comm_ig2_b=2.177,
              fpga_bias_id1=18.008, fpga_bias_id2=12.998)
    
    print('\n\n\nHK\n\n\n')
    
    for CCD in CCDs:
        CCDk = 'CCD%i' % CCD
        for Q in Quads:
            
            for mkey in mkeys:
                
                ftrans_ELVIS = ELVIS_translators[mkey]  
                
                HKkey = ftrans_ELVIS(CCD,Q,mode='HK')
                
                EVHK = HK[HKkey]
                
                countsHK = f_V2counts_EHK(EVHK,HKkey,version='7.5.X')
                
                ftrans_CAL = V_translators[mkey]
                Vkey = ftrans_CAL(CCD,Q)
                
                if mkey in ['IDL','IDH']:
                    _mkey = 'ID'
                else:
                    _mkey = mkey
                
                keytuple = (CCDk,Q,_mkey)
                
                trueVHK = f_counts2V_CAL(countsHK, VCALdict, keytuple, mode='ADC')
                
                print('%s=%.1f     %s=%.2f' % (HKkey, EVHK, Vkey, trueVHK))
                
    

if __name__ == '__main__':
    """ """
    
    calpath = 'FPA_CAL/Boxed_FM_Units/CCD_CAL'
    calfile = os.path.join(calpath,'CCD_CAL_CONVERSIONS_ALL_BLOCKS.xlsx')
    ROE=3
    
    proto(calfile,ROE)
    
    