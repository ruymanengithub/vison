#! /unsafe/raf/SOFTWARE/anaconda/envs/VISSIM/bin/python
# -*- coding: utf-8 -*-
"""

House-Keeping inspection and handling tools.

:History:
Created on Thu Mar 10 12:11:58 2016

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import os
from glob import glob
import datetime
import sys
from collections import OrderedDict
import copy

from astropy.io import ascii
from astropy.table import Table, Column, vstack
import string as st
import numpy as np

from matplotlib import pyplot as plt
import matplotlib
matplotlib.rcParams['font.size'] = 17
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('axes', linewidth=1.1)
matplotlib.rcParams['legend.fontsize'] = 12
matplotlib.rcParams['legend.handlelength'] = 3
matplotlib.rcParams['xtick.major.size'] = 5
matplotlib.rcParams['ytick.major.size'] = 5
matplotlib.rcParams['image.interpolation'] = 'none'
import pylab
#import matplotlib.cm as cm
#from mpl_toolkits.mplot3d import Axes3D

from vison.support import context, utils
#from vison.pipe import lib as pilib


# END IMPORT


allHK_keys = {
    '5.7.02': [
        'HK_OD_Top_CCD1',
        'HK_OD_Bottom_CCD1',
        'HK_OD_Top_CCD2',
        'HK_OD_Bottom_CCD2',
        'HK_OD_Top_CCD3',
        'HK_OD_Bottom_CCD3',
        'HK_RD_Top_CCD1',
        'HK_RD_Bottom_CCD1',
        'HK_RD_Top_CCD2',
        'HK_RD_Bottom_CCD2',
        'HK_RD_Top_CCD3',
        'HK_RD_Bottom_CCD3',
        'HK_temp_top_CCD1',
        'HK_temp_bottom_CCD1',
        'HK_temp_top_CCD2',
        'HK_temp_bottom_CCD2',
        'HK_temp_top_CCD3',
        'HK_temp_bottom_CCD3',
        'HK_IG1_top',
        'HK_IG1_bot',
        'HK_IG2_top',
        'HK_IG2_bot',
        'HK_IDH',
        'HK_IDL',
        'HK_DD_bias',
        'HK_OG_bias',
        'HK_1.5V_ROE',
        'HK_VCCD_ROE',
        'HK_5VA_pos_ROE',
        'HK_5V_ref_ROE',
        'HK_10VA_ROE',
        'HK_5.2V_neg_ROE',
        'HK_3V_neg_ROE',
        'HK_VRclk_ROE',
        'HK_VRClk_Lo_ROE',
        'HK_3.3V_DIG_RPSU',
        'HK_I3.3V_DIG_RPSU',
        'HK_1.5V_DIG_RPSU',
        'HK_I1.5V_DIG_RPSU',
        'HK_28V_Pri_RPSU',
        'HK_I28V_RPSU',
        'HK_VAN_pos_RPSU',
        'HK_I+VAN_RPSU',
        'HK_VAN_neg_RPSU',
        'HK_I-VAN_RPSU',
        'HK_VCLK_RPSU',
        'HK_IVCLK_RPSU',
        'HK_VCCD_RPSU',
        'HK_IVCCD_RPSU',
        'HK_Temp1_RPSU',
        'HK_3.3V_ROE',
        'HK_Temp2_RPSU',
        'HK_VClk_ROE'],
    '5.7.04': [
        'HK_OD_Top_CCD1',
        'HK_OD_Bottom_CCD1',
        'HK_OD_Top_CCD2',
        'HK_OD_Bottom_CCD2',
        'HK_OD_Top_CCD3',
        'HK_OD_Bottom_CCD3',
        'HK_RD_Top_CCD1',
        'HK_RD_Bottom_CCD1',
        'HK_RD_Top_CCD2',
        'HK_RD_Bottom_CCD2',
        'HK_RD_Top_CCD3',
        'HK_RD_Bottom_CCD3',
        'HK_temp_top_CCD1',
        'HK_temp_bottom_CCD1',
        'HK_temp_top_CCD2',
        'HK_temp_bottom_CCD2',
        'HK_temp_top_CCD3',
        'HK_temp_bottom_CCD3',
        'HK_IG1_top',
        'HK_IG1_bot',
        'HK_IG2_top',
        'HK_IG2_bot',
        'HK_IDH',
        'HK_IDL',
        'HK_DD_bias',
        'HK_OG_bias',
        'HK_1.5V_ROE',
        'HK_VCCD_ROE',
        'HK_5VA_pos_ROE',
        'HK_5V_ref_ROE',
        'HK_10VA_ROE',
        'HK_5.2V_neg_ROE',
        'HK_3V_neg_ROE',
        'HK_VRclk_ROE',
        'HK_VRClk_Lo_ROE',
        'HK_3.3V_DIG_RPSU',
        'HK_I3.3V_DIG_RPSU',
        'HK_1.5V_DIG_RPSU',
        'HK_I1.5V_DIG_RPSU',
        'HK_28V_Pri_RPSU',
        'HK_I28V_RPSU',
        'HK_VAN_pos_RPSU',
        'HK_I+VAN_RPSU',
        'HK_VAN_neg_RPSU',
        'HK_I-VAN_RPSU',
        'HK_VCLK_RPSU',
        'HK_IVCLK_RPSU',
        'HK_VCCD_RPSU',
        'HK_IVCCD_RPSU',
        'HK_Temp1_RPSU',
        'HK_3.3V_ROE',
        'HK_Temp2_RPSU',
        'HK_VClk_ROE'],
    '5.7.07': [
        'HK_OD_Top_CCD1',
        'HK_OD_Bottom_CCD1',
        'HK_OD_Top_CCD2',
        'HK_OD_Bottom_CCD2',
        'HK_OD_Top_CCD3',
        'HK_OD_Bottom_CCD3',
        'HK_RD_Top_CCD1',
        'HK_RD_Bottom_CCD1',
        'HK_RD_Top_CCD2',
        'HK_RD_Bottom_CCD2',
        'HK_RD_Top_CCD3',
        'HK_RD_Bottom_CCD3',
        'HK_temp_top_CCD1',
        'HK_temp_bottom_CCD1',
        'HK_temp_top_CCD2',
        'HK_temp_bottom_CCD2',
        'HK_temp_top_CCD3',
        'HK_temp_bottom_CCD3',
        'HK_IG1_top',
        'HK_IG1_bot',
        'HK_IG2_top',
        'HK_IG2_bot',
        'HK_IDH',
        'HK_IDL',
        'HK_DD_bias',
        'HK_OG_bias',
        'HK_1.5V_ROE',
        'HK_VCCD_ROE',
        'HK_5VA_pos_ROE',
        'HK_5V_ref_ROE',
        'HK_10VA_ROE',
        'HK_5.2V_neg_ROE',
        'HK_3V_neg_ROE',
        'HK_VRclk_ROE',
        'HK_VRClk_Lo_ROE',
        'HK_3.3V_DIG_RPSU',
        'HK_I3.3V_DIG_RPSU',
        'HK_1.5V_DIG_RPSU',
        'HK_I1.5V_DIG_RPSU',
        'HK_28V_Pri_RPSU',
        'HK_I28V_RPSU',
        'HK_VAN_pos_RPSU',
        'HK_I+VAN_RPSU',
        'HK_VAN_neg_RPSU',
        'HK_I-VAN_RPSU',
        'HK_VCLK_RPSU',
        'HK_IVCLK_RPSU',
        'HK_VCCD_RPSU',
        'HK_IVCCD_RPSU',
        'HK_Temp1_RPSU',
        'HK_3.3V_ROE',
        'HK_Temp2_RPSU',
        'HK_VClk_ROE',
        'Video_TOP',
        'Video_BOT']}

allHK_keys['5.7.06'] = allHK_keys['5.7.04']
allHK_keys['5.7.08'] = allHK_keys['5.7.07']
allHK_keys['5.7.09'] = allHK_keys['5.7.08']

#allHK_keys['5.8.X'] = allHK_keys['5.7.09']
#allHK_keys['5.8.X'] = ['HK_timestamp'] + allHK_keys['5.8.X']

allHK_keys['6.0.0'] = ['TimeStamp', 'HK_OD_Top_CCD1', 'HK_OD_Bottom_CCD1',
                       'HK_OD_Top_CCD2', 'HK_OD_Bottom_CCD2', 'HK_OD_Top_CCD3', 'HK_OD_Bottom_CCD3',
                       'HK_IG1_Top_CCD1', 'HK_IG1_Bottom_CCD1', 'HK_IG1_Top_CCD2', 'HK_IG1_Bottom_CCD2',
                       'HK_IG1_Top_CCD3', 'HK_IG1_Bottom_CCD3', 'HK_temp_top_CCD1', 'HK_temp_bottom_CCD1',
                       'HK_temp_top_CCD2', 'HK_temp_bottom_CCD2', 'HK_temp_top_CCD3',
                       'HK_temp_bottom_CCD3', 'HK_RD_top', 'HK_RD_bot', 'HK_IG2_top', 'HK_IG2_bot',
                       'HK_IDH', 'HK_IDL', 'HK_DD_bias', 'HK_OG_bias', 'HK_1.5V_ROE', 'HK_VCCD_ROE',
                       'HK_5VA_pos_ROE', 'HK_5V_ref_ROE', 'HK_10VA_ROE', 'HK_5.2V_neg_ROE', 'HK_3V_neg_ROE',
                       'HK_VRclk_ROE', 'HK_VRClk_Lo_ROE', 'HK_3.3V_DIG_RPSU', 'HK_I3.3V_DIG_RPSU',
                       'HK_1.5V_DIG_RPSU', 'HK_I1.5V_DIG_RPSU', 'HK_28V_Pri_RPSU', 'HK_I28V_RPSU',
                       'HK_VAN_pos_RPSU', 'HK_I+VAN_RPSU', 'HK_VAN_neg_RPSU', 'HK_I-VAN_RPSU',
                       'HK_VCLK_RPSU', 'HK_IVCLK_RPSU', 'HK_VCCD_RPSU', 'HK_IVCCD_RPSU', 'HK_Temp1_RPSU',
                       'HK_Temp2_RPSU', 'HK_Video_TOP', 'HK_Video_BOT', 'HK_FPGA_TOP', 'HK_FPGA_BOT',
                       'HK_ID1', 'HK_ID2', 'HK_Viclk_ROE']

allHK_keys['6.1.0'] = allHK_keys['6.0.0']

allHK_keys['6.3.0'] = [
    'TimeStamp',
    'CCD1_OD_T',
    'CCD2_OD_T',
    'CCD3_OD_T',
    'COMM_RD_T',
    'CCD2_IG1_T',
    'CCD3_IG1_T',
    'CCD1_TEMP_T',
    'CCD2_TEMP_T',
    'CCD3_TEMP_T',
    'CCD1_IG1_T',
    'COMM_IG2_T',
    'VID_PCB_TEMP_T',
    'FPGA_PCB_TEMP_T',
    'CCD1_OD_B',
    'CCD2_OD_B',
    'CCD3_OD_B',
    'COMM_RD_B',
    'CCD2_IG1_B',
    'CCD3_IG1_B',
    'CCD1_TEMP_B',
    'CCD2_TEMP_B',
    'CCD3_TEMP_B',
    'CCD1_IG1_B',
    'COMM_IG2_B',
    'VID_PCB_TEMP_B',
    'FPGA_PCB_TEMP_B',
    'FPGA_BIAS_DD',
    'FPGA_BIAS_OG',
    'FPGA_BIAS_ID1',
    'FPGA_BIAS_ID2',
    'FPGA_BIAS_ID_T',
    'FPGA_BIAS_ID_B',
    'FPGA_VRCLK_V',
    'FPGA_10VA_P_V',
    'FPGA_VICLK_V',
    'FPGA_5VA_P_V',
    'FPGA_5VREF_V',
    'FPGA_VCCD_V',
    'FPGA_1V5VD_P_V',
    'FPGA_3V2_N_V',
    'FPGA_5VA_N_V',
    'RPSU_VCCD_V',
    'RPSU_VCLK_V',
    'RPSU_VAN_P_V',
    'RPSU_VAN_N_V',
    'RPSU_3V3VD_V',
    'RPSU_1V5VD_V',
    'RPSU_28V_PRI_V',
    'RPSU_TEMP1',
    'RPSU_VCCD_I',
    'RPSU_VCLK_I',
    'RPSU_VAN_P_I',
    'RPSU_VAN_N_I',
    'RPSU_3V3VD_I',
    'RPSU_1V5VD_I',
    'RPSU_28V_PRI_I',
    'RPSU_TEMP_2',
    'stTMRErrFlg',
    'hkTMRErrFlg',
    'spwTmTOFlg',
    'CDPUClkSt',
    'fpgaSpwErr',
    'V3v3ProtCnt',
    'V5ProtCnt',
    'VccdErrFlg',
    'VclkErrFlg',
    'VanErrFlg',
    'stRamErr',
    'hkRamErr',
    'ADC_BSY_ERR_CNT',
    'SPW_STATUS_REG']


allHK_keys['6.5.X'] = ['TimeStamp', 'CCD3_OD_T', 'CCD2_OD_T', 'CCD1_OD_T',
                       'COMM_RD_T', 'CCD2_IG1_T', 'CCD1_IG1_T',
                       'CCD3_TEMP_T', 'CCD2_TEMP_T', 'CCD1_TEMP_T', 'CCD3_IG1_T', 'COMM_IG2_T',
                       'VID_PCB_TEMP_T', 'FPGA_PCB_TEMP_T',
                       'CCD1_OD_B', 'CCD2_OD_B', 'CCD3_OD_B',
                       'COMM_RD_B', 'CCD2_IG1_B', 'CCD3_IG1_B',
                       'CCD1_TEMP_B', 'CCD2_TEMP_B', 'CCD3_TEMP_B',
                       'CCD1_IG1_B', 'COMM_IG2_B',
                       'VID_PCB_TEMP_B', 'FPGA_PCB_TEMP_B',
                       'FPGA_BIAS_DD', 'FPGA_BIAS_OG', 'FPGA_BIAS_ID1', 'FPGA_BIAS_ID2',
                       'FPGA_BIAS_ID_T', 'FPGA_BIAS_ID_B',
                       'FPGA_VRCLK_V', 'FPGA_10VA_P_V', 'FPGA_VICLK_V', 'FPGA_5VA_P_V',
                       'FPGA_5VREF_V', 'FPGA_VCCD_V', 'FPGA_1V5VD_P_V', 'FPGA_3V2_N_V',
                       'FPGA_5VA_N_V',
                       'RPSU_VCCD_V', 'RPSU_VCLK_V', 'RPSU_VAN_P_V', 'RPSU_VAN_N_V', 'RPSU_3V3VD_V',
                       'RPSU_1V5VD_V', 'RPSU_28V_PRI_V', 'RPSU_TEMP1', 'RPSU_VCCD_I', 'RPSU_VCLK_I',
                       'RPSU_VAN_P_I', 'RPSU_VAN_N_I', 'RPSU_3V3VD_I', 'RPSU_1V5VD_I',
                       'RPSU_28V_PRI_I', 'RPSU_TEMP_2',
                       'ProtOvRideFlg', 'hkInvalidFlg', 'SPI_Inh_n', '3v3ProtErr', '5vProtErr',
                       'V3v3ProtCnt', 'V5ProtCnt', 'VccdErrFlg', 'VclkErrFlg', 'VanErrFlg',
                       'hkTmErr', 'spwTmTOFlg', 'CDPUClkSt', 'fpgaSpwErr', 'hkRamErr',
                       'ADC_BSY_ERR_CNT', 'SPW_STATUS_REG', 'Out_of_range']

allHK_keys['7.2.X'] = ['TimeStamp', 'CCD3_OD_T', 'CCD2_OD_T', 'CCD1_OD_T', 'COMM_RD_T',
                       'CCD2_IG1_T', 'CCD1_IG1_T', 'CCD3_TEMP_T', 'CCD2_TEMP_T', 'CCD1_TEMP_T', 'CCD3_IG1_T',
                       'COMM_IG2_T', 'VID_PCB_TEMP_T', 'FPGA_PCB_TEMP_T', 'CCD1_OD_B', 'CCD2_OD_B',
                       'CCD3_OD_B', 'COMM_RD_B', 'CCD2_IG1_B', 'CCD3_IG1_B', 'CCD1_TEMP_B',
                       'CCD2_TEMP_B', 'CCD3_TEMP_B', 'CCD1_IG1_B', 'COMM_IG2_B', 'VID_PCB_TEMP_B',
                       'FPGA_PCB_TEMP_B', 'FPGA_BIAS_DD', 'FPGA_BIAS_OG', 'FPGA_BIAS_ID1', 'FPGA_BIAS_ID2',
                       'FPGA_BIAS_ID_T', 'FPGA_BIAS_ID_B', 'FPGA_VRCLK_V', 'FPGA_10VA_P_V', 'FPGA_VICLK_V',
                       'FPGA_5VA_P_V', 'FPGA_5VREF_V', 'FPGA_VCCD_V', 'FPGA_1V5VD_P_V', 'FPGA_3V2_N_V',
                       'FPGA_5VA_N_V', 'RPSU_VCCD_V', 'RPSU_VCLK_V', 'RPSU_VAN_P_V', 'RPSU_VAN_N_V',
                       'RPSU_3V3VD_V', 'RPSU_1V5VD_V', 'RPSU_28V_PRI_V', 'RPSU_TEMP1', 'RPSU_VCCD_I',
                       'RPSU_VCLK_I', 'RPSU_VAN_P_I', 'RPSU_VAN_N_I', 'RPSU_3V3VD_I', 'RPSU_1V5VD_I',
                       'RPSU_28V_PRI_I', 'RPSU_TEMP_2', 'ProtOvRideFlg', 'hkInvalidFlg', 'SPI_Inh_n',
                       '3v3ProtErr', '5vProtErr', 'V3v3ProtCnt', 'V5ProtCnt', 'VccdErrFlg', 'VclkErrFlg',
                       'VanErrFlg', 'dacExeErr', 'hkTmErr', 'spwTmTOFlg', 'spwTmRspErrFlg', 'spwTcPktErrFlg',
                       'adcBsyErrFlg',
                       'CDPUClkSt', 'CDPUClkLost', 'hkRamErr', 'ADC_BSY_ERR_CNT', 'SPW_STATUS_REG',
                       'Out_of_range']

allHK_keys['7.5.X'] = ['TimeStamp', 'CCD3_OD_T', 'CCD2_OD_T', 'CCD1_OD_T', 'COMM_RD_T',
                       'CCD2_IG1_T', 'CCD1_IG1_T', 'CCD3_TEMP_T', 'CCD2_TEMP_T', 'CCD1_TEMP_T', 'CCD3_IG1_T',
                       'COMM_IG2_T', 'VID_PCB_TEMP_T', 'FPGA_PCB_TEMP_T', 'CCD1_OD_B', 'CCD2_OD_B',
                       'CCD3_OD_B', 'COMM_RD_B', 'CCD2_IG1_B', 'CCD3_IG1_B', 'CCD1_TEMP_B',
                       'CCD2_TEMP_B', 'CCD3_TEMP_B', 'CCD1_IG1_B', 'COMM_IG2_B', 'VID_PCB_TEMP_B',
                       'FPGA_PCB_TEMP_B', 'FPGA_BIAS_DD', 'FPGA_BIAS_OG', 'FPGA_BIAS_ID1', 'FPGA_BIAS_ID2',
                       'FPGA_BIAS_ID_T', 'FPGA_BIAS_ID_B', 'FPGA_VRCLK_V', 'FPGA_10VA_P_V', 'FPGA_VICLK_V',
                       'FPGA_5VA_P_V', 'FPGA_5VREF_V', 'FPGA_VCCD_V', 'FPGA_1V5VD_P_V', 'FPGA_3V2_N_V',
                       'FPGA_5VA_N_V', 'RPSU_VCCD_V', 'RPSU_VCLK_V', 'RPSU_VAN_P_V', 'RPSU_VAN_N_V',
                       'RPSU_3V3VD_V', 'RPSU_1V5VD_V', 'RPSU_28V_PRI_V', 'RPSU_TEMP1', 'RPSU_VCCD_I',
                       'RPSU_VCLK_I', 'RPSU_VAN_P_I', 'RPSU_VAN_N_I', 'RPSU_3V3VD_I', 'RPSU_1V5VD_I',
                       'RPSU_28V_PRI_I', 'RPSU_TEMP_2', 'ProtOvRideFlg', 'hkInvalidFlg', 'SPI_Inh_n',
                       '3v3ProtErr', '5vProtErr', 'V3v3ProtCnt', 'V5ProtCnt', 'VccdErrFlg', 'VclkErrFlg',
                       'VanErrFlg', 'dacExeErr', 'hkTmErr', 'spwTmTOFlg', 'spwTmRspErrFlg', 'spwTcPktErrFlg',
                       'adcBsyErrFlg',
                       'CDPUClkSw', 'hkRamErr', 'ADC_BSY_ERR_CNT', 'pixelToErr', 'spwStatus',
                       'spwTmPtypeFlag', 'TwoCmdsErr', 'SPW_STATUS_REG', 'Out_of_range']

# HKlims = dict(Performance=dict(key=[Relative/Absolute/Identity,min,max]))

HKlims = {}
HKlims['6.3.0'] = dict(P={'CCD1_OD_T': ['R', -0.2, +0.2], 'CCD2_OD_T': ['R', -0.2, +0.2], 'CCD3_OD_T': ['R', -0.2, +0.2],
                          'CCD1_OD_B': ['R', -0.2, +0.2], 'CCD2_OD_B': ['R', -0.2, +0.2], 'CCD3_OD_B': ['R', -0.2, +0.2],
                          'COMM_RD_T': ['R', -0.2, +0.2], 'COMM_RD_B': ['R', -0.2, +0.2],
                          'CCD1_IG1_T': ['R', -0.2, +0.2], 'CCD2_IG1_T': ['R', -0.2, +0.2], 'CCD3_IG1_T': ['R', -0.2, +0.2],
                          'CCD2_IG1_B': ['R', -0.2, +0.2], 'CCD3_IG1_B': ['R', -0.2, +0.2], 'CCD1_IG1_B': ['R', -0.2, +0.2],
                          'COMM_IG2_T': ['R', -0.2, +0.2], 'COMM_IG2_B': ['R', -0.2, +0.2],
                          'CCD1_TEMP_T': ['A', -121., -119.], 'CCD2_TEMP_T': ['A', -121., -119.], 'CCD3_TEMP_T': ['A', -121., -119.],
                          'CCD1_TEMP_B': ['A', -121., -119.], 'CCD2_TEMP_B': ['A', -121., -119.], 'CCD3_TEMP_B': ['A', -121., -119.],
                          'VID_PCB_TEMP_T': ['A', -10., +40.], 'FPGA_PCB_TEMP_T': ['A', -10., +40.],
                          'VID_PCB_TEMP_B': ['A', -10., +40.], 'FPGA_PCB_TEMP_B': ['A', -10., +40.],
                          'FPGA_BIAS_DD': ['A', 21.8, 22.2], 'FPGA_BIAS_OG': ['A', 1.8, 2.2],
                          'FPGA_BIAS_ID1': ['R', -0.2, +0.2], 'FPGA_BIAS_ID2': ['R', -0.2, +0.2],
                          'FPGA_BIAS_ID_T': ['A', -0.3, +25.], 'FPGA_BIAS_ID_B': ['A', -0.3, +25.],
                          'FPGA_VRCLK_V': ['A', 9., 11.], 'FPGA_10VA_P_V': ['A', 9., 11.],
                          'FPGA_VICLK_V': ['A', 7.2, 8.8], 'FPGA_5VA_P_V': ['A', 4.75, 5.25], 'FPGA_5VREF_V': ['A', 4.9, 5.1],
                          'FPGA_VCCD_V': ['A', 29.5, 30.5], 'FPGA_1V5VD_P_V': ['A', 1.45, 1.55],
                          'FPGA_3V2_N_V': ['A', -3.7, -2.7], 'FPGA_5VA_N_V': ['A', -5.9, -4.9],
                          'RPSU_VCCD_V': ['A', 31.6, 36.], 'RPSU_VCLK_V': ['A', 11.5, 14.], 'RPSU_VAN_P_V': ['A', 6.0, 8.2],
                          'RPSU_VAN_N_V': ['A', -8.2, -6.], 'RPSU_3V3VD_V': ['A', 3.1, 3.6], 'RPSU_1V5VD_V': ['A', 1.435, 1.575],
                          'RPSU_28V_PRI_V': ['A', 26., 29.], 'RPSU_TEMP1': ['A', -10., 40.],
                          'RPSU_VCCD_I': ['A', 0., 0.1], 'RPSU_VCLK_I': ['A', 0., 0.1], 'RPSU_VAN_P_I': ['A', 0., 0.6],
                          'RPSU_VAN_N_I': ['A', 0.0, 0.6], 'RPSU_3V3VD_I': ['A', 0., 0.150], 'RPSU_1V5VD_I': ['A', 0.0, 0.150],
                          'RPSU_28V_PRI_I': ['A', 0.35, 0.5], 'RPSU_TEMP_2': ['A', -10., 40.],
                          'stTMRErrFlg': ['I', 0], 'hkTMRErrFlg': ['I', 0],
                          'spwTmTOFlg': ['I', 0], 'CDPUClkSt': ['I', 0], 'fpgaSpwErr': ['I', 0], 'V3v3ProtCnt': ['I', 0], 'V5ProtCnt': ['I', 0],
                          'VccdErrFlg': ['I', 0], 'VclkErrFlg': ['I', 0], 'VanErrFlg': ['I', 0], 'stRamErr': ['I', 0], 'hkRamErr': ['I', 0],
                          'ADC_BSY_ERR_CNT': ['I', 0],
                          'SPW_STATUS_REG': ['I', '2FA0']},
                       S={'CCD1_OD_T': ['A', -0.3, +30.], 'CCD2_OD_T': ['A', -0.3, +30.], 'CCD3_OD_T': ['A', -0.3, +30.],
                          'CCD1_OD_B': ['A', -0.3, +30.], 'CCD2_OD_B': ['A', -0.3, +30.], 'CCD3_OD_B': ['A', -0.3, +30.],
                          'COMM_RD_T': ['A', -0.3, +25.], 'COMM_RD_B': ['A', -0.3, +25.],
                          'CCD1_IG1_T': ['A', -20., +20.], 'CCD2_IG1_T': ['A', -20., +20.], 'CCD3_IG1_T': ['A', -20., +20.],
                          'CCD2_IG1_B': ['A', -20., +20.], 'CCD3_IG1_B': ['A', -20., +20.], 'CCD1_IG1_B': ['A', -20., +20.],
                          'COMM_IG2_T': ['A', -20., +20.], 'COMM_IG2_B': ['A', -20., +20.],
                          'CCD1_TEMP_T': ['A', -133., 35.], 'CCD2_TEMP_T': ['A', -133., 35.], 'CCD3_TEMP_T': ['A', -133., 35.],
                          'CCD1_TEMP_B': ['A', -133., 35.], 'CCD2_TEMP_B': ['A', -133., 35.], 'CCD3_TEMP_B': ['A', -133., 35.],
                          'VID_PCB_TEMP_T': ['A', -10., +40.], 'FPGA_PCB_TEMP_T': ['A', -10., +40.],
                          'VID_PCB_TEMP_B': ['A', -10., +40.], 'FPGA_PCB_TEMP_B': ['A', -10., +40.],
                          'FPGA_BIAS_DD': ['A', -0.3, 30.], 'FPGA_BIAS_OG': ['A', -20., 20.],
                          'FPGA_BIAS_ID1': ['A', -0.3, +25.], 'FPGA_BIAS_ID2': ['A', -0.3, +25.],
                          'FPGA_BIAS_ID_T': ['A', -0.3, +25.], 'FPGA_BIAS_ID_B': ['A', -0.3, +25.],
                          'FPGA_VRCLK_V': ['A', 9., 11.], 'FPGA_10VA_P_V': ['A', 9., 11.],
                          'FPGA_VICLK_V': ['A', 7.2, 8.8], 'FPGA_5VA_P_V': ['A', 4.75, 5.25], 'FPGA_5VREF_V': ['A', 4.9, 5.1],
                          'FPGA_VCCD_V': ['A', 29.5, 30.5], 'FPGA_1V5VD_P_V': ['A', 1.45, 1.55],
                          'FPGA_3V2_N_V': ['A', -3.7, -2.7], 'FPGA_5VA_N_V': ['A', -5.9, -4.9],
                          'RPSU_VCCD_V': ['A', 31.6, 36.], 'RPSU_VCLK_V': ['A', 11.5, 14.], 'RPSU_VAN_P_V': ['A', 6.0, 8.2],
                          'RPSU_VAN_N_V': ['A', -8.2, -6.], 'RPSU_3V3VD_V': ['A', 3.1, 3.6], 'RPSU_1V5VD_V': ['A', 1.435, 1.575],
                          'RPSU_28V_PRI_V': ['A', 26., 29.], 'RPSU_TEMP1': ['A', -10., 40.],
                          'RPSU_VCCD_I': ['A', 0., 0.1], 'RPSU_VCLK_I': ['A', 0., 0.1], 'RPSU_VAN_P_I': ['A', 0., 0.6],
                          'RPSU_VAN_N_I': ['A', 0.0, 0.6], 'RPSU_3V3VD_I': ['A', 0., 0.150], 'RPSU_1V5VD_I': ['A', 0.0, 0.150],
                          'RPSU_28V_PRI_I': ['A', 0.35, 0.5], 'RPSU_TEMP_2': ['A', -10., 40.],
                          'stTMRErrFlg': ['I', 0], 'hkTMRErrFlg': ['I', 0],
                          'spwTmTOFlg': ['I', 0], 'CDPUClkSt': ['I', 0], 'fpgaSpwErr': ['I', 0], 'V3v3ProtCnt': ['I', 0], 'V5ProtCnt': ['I', 0],
                          'VccdErrFlg': ['I', 0], 'VclkErrFlg': ['I', 0], 'VanErrFlg': ['I', 0], 'stRamErr': ['I', 0], 'hkRamErr': ['I', 0],
                          'ADC_BSY_ERR_CNT': ['I', 0],
                          'SPW_STATUS_REG': ['I', '2FA0']})

HKlims['6.5.X'] = dict()
HKlims['6.5.X']['P'] = HKlims['6.3.0']['P']
ignore = list(map(HKlims['6.5.X']['P'].pop,
             ['stTMRErrFlg',
              'hkTMRErrFlg',
              'spwTmTOFlg',
              'CDPUClkSt',
              'fpgaSpwErr',
              'V3v3ProtCnt',
              'V5ProtCnt',
              'VccdErrFlg',
              'VclkErrFlg',
              'VanErrFlg',
              'stRamErr',
              'hkRamErr',
              'ADC_BSY_ERR_CNT',
              'SPW_STATUS_REG']))

HKlims['6.5.X']['P'].update(
    {
        'ProtOvRideFlg': [
            'I', 0], 'hkInvalidFlg': [
                'I', 0], 'SPI_Inh_n': [
                    'I', 0], '3v3ProtErr': [
                        'I', 0], '5vProtErr': [
                            'I', 0], 'V3v3ProtCnt': [
                                'I', 0], 'V5ProtCnt': [
                                    'I', 0], 'VccdErrFlg': [
                                        'I', 0], 'VclkErrFlg': [
                                            'I', 0], 'VanErrFlg': [
                                                'I', 0], 'hkTmErr': [
                                                    'I', 0], 'spwTmTOFlg': [
                                                        'I', 0], 'CDPUClkSt': [
                                                            'I', 0], 'fpgaSpwErr': [
                                                                'I', 0], 'hkRamErr': [
                                                                    'I', 0], 'ADC_BSY_ERR_CNT': [
                                                                        'I', 0], 'SPW_STATUS_REG': [
                                                                            'I', 12192], 'Out_of_range': [
                                                                                'I', 0]})

HKlims['6.5.X']['S'] = HKlims['6.3.0']['S']
ignore = list(map(HKlims['6.5.X']['S'].pop,
             ['stTMRErrFlg',
              'hkTMRErrFlg',
              'spwTmTOFlg',
              'CDPUClkSt',
              'fpgaSpwErr',
              'V3v3ProtCnt',
              'V5ProtCnt',
              'VccdErrFlg',
              'VclkErrFlg',
              'VanErrFlg',
              'stRamErr',
              'hkRamErr',
              'ADC_BSY_ERR_CNT',
              'SPW_STATUS_REG']))

HKlims['6.5.X']['S'].update(
    {
        'ProtOvRideFlg': [
            'I', 0], 'hkInvalidFlg': [
                'I', 0], 'SPI_Inh_n': [
                    'I', 0], '3v3ProtErr': [
                        'I', 0], '5vProtErr': [
                            'I', 0], 'V3v3ProtCnt': [
                                'I', 0], 'V5ProtCnt': [
                                    'I', 0], 'VccdErrFlg': [
                                        'I', 0], 'VclkErrFlg': [
                                            'I', 0], 'VanErrFlg': [
                                                'I', 0], 'hkTmErr': [
                                                    'I', 0], 'spwTmTOFlg': [
                                                        'I', 0], 'CDPUClkSt': [
                                                            'I', 0], 'fpgaSpwErr': [
                                                                'I', 0], 'hkRamErr': [
                                                                    'I', 0], 'ADC_BSY_ERR_CNT': [
                                                                        'I', 0], 'SPW_STATUS_REG': [
                                                                            'I', 12192], 'Out_of_range': [
                                                                                'I', 0]})

HKlims['7.2.X'] = dict()
HKlims['7.2.X']['P'] = HKlims['6.5.X']['P'].copy()
ignore = list(map(HKlims['7.2.X']['P'].pop, ['fpgaSpwErr']))
HKlims['7.2.X']['P'].update({'dacExeErr': ['I', 0], 'spwTmRspErrFlg': ['I', 0], 'spwTcPktErrFlg': [
                            'I', 0], 'adcBsyErrFlg': ['I', 0], 'CDPUClkLost': ['I', 0]})

HKlims['7.2.X']['S'] = HKlims['6.5.X']['S'].copy()
ignore = list(map(HKlims['7.2.X']['S'].pop, ['fpgaSpwErr']))
HKlims['7.2.X']['S'].update({'dacExeErr': ['I', 0], 'spwTmRspErrFlg': ['I', 0], 'spwTcPktErrFlg': [
                            'I', 0], 'adcBsyErrFlg': ['I', 0], 'CDPUClkLost': ['I', 0]})

HKlims['7.5.X'] = dict()

for q in ['P', 'S']:
    HKlims['7.5.X'][q] = HKlims['7.2.X'][q].copy()
    ignore = list(map(HKlims['7.2.X'][q].pop, ['CDPUClkLost', 'CDPUClkSt']))
    HKlims['7.5.X'][q].update(
        dict(
            CDPUClkSw=[
                'I', 0], pixelToErr=[
                'I', 0], spwStatus=[
                    'I', 0], spwTmPtypeFlag=[
                        'I', 0], TwoCmdsErr=[
                            'I', 0]))


HKcorr = {}
HKcorr['6.3.0'] = {
    'CCD1_OD_T': 'OD_1_T',
    'CCD2_OD_T': 'OD_2_T',
    'CCD3_OD_T': 'OD_3_T',
    'CCD1_OD_B': 'OD_1_B',
    'CCD2_OD_B': 'OD_2_B',
    'CCD3_OD_B': 'OD_3_B',
    'COMM_RD_T': 'RD_T',
    'COMM_RD_B': 'RD_B',
    'CCD1_IG1_T': 'IG1_1_T',
    'CCD2_IG1_T': 'IG1_2_T',
    'CCD3_IG1_T': 'IG1_3_T',
    'CCD2_IG1_B': 'IG1_2_T',
    'CCD3_IG1_B': 'IG1_3_T',
    'CCD1_IG1_B': 'IG1_1_B',
    'COMM_IG2_T': 'IG2_T',
    'COMM_IG2_B': 'IG2_B',
    'CCD1_TEMP_T': None,
    'CCD2_TEMP_T': None,
    'CCD3_TEMP_T': None,
    'CCD1_TEMP_B': None,
    'CCD2_TEMP_B': None,
    'CCD3_TEMP_B': None,
    'VID_PCB_TEMP_T': None,
    'FPGA_PCB_TEMP_T': None,
    'VID_PCB_TEMP_B': None,
    'FPGA_PCB_TEMP_B': None,
    'FPGA_BIAS_DD': None,
    'FPGA_BIAS_OG': None,
    'FPGA_BIAS_ID1': 'IDH',
    'FPGA_BIAS_ID2': 'IDL',
    'FPGA_BIAS_ID_T': None,
    'FPGA_BIAS_ID_B': None,
    'FPGA_VRCLK_V': None,
    'FPGA_10VA_P_V': None,
    'FPGA_VICLK_V': None,
    'FPGA_5VA_P_V': None,
    'FPGA_5VREF_V': None,
    'FPGA_VCCD_V': None,
    'FPGA_1V5VD_P_V': None,
    'FPGA_3V2_N_V': None,
    'FPGA_5VA_N_V': None,
    'RPSU_VCCD_V': None,
    'RPSU_VCLK_V': None,
    'RPSU_VAN_P_V': None,
    'RPSU_VAN_N_V': None,
    'RPSU_3V3VD_V': None,
    'RPSU_1V5VD_V': None,
    'RPSU_28V_PRI_V': None,
    'RPSU_TEMP1': None,
    'RPSU_VCCD_I': None,
    'RPSU_VCLK_I': None,
    'RPSU_VAN_P_I': None,
    'RPSU_VAN_N_I': None,
    'RPSU_3V3VD_I': None,
    'RPSU_1V5VD_I': None,
    'RPSU_28V_PRI_I': None,
    'RPSU_TEMP_2': None,
    'ProtOvRideFlg': None,
    'hkInvalidFlg': None,
    'SPI_Inh_n': None,
    '3v3ProtErr': None,
    '5vProtErr': None,
    'V3v3ProtCnt': None,
    'V5ProtCnt': None,
    'VccdErrFlg': None,
    'VclkErrFlg': None,
    'VanErrFlg': None,
    'hkTmErr': None,
    'spwTmTOFlg': None,
    'CDPUClkSt': None,
    'fpgaSpwErr': None,
    'hkRamErr': None,
    'ADC_BSY_ERR_CNT': None,
    'SPW_STATUS_REG': None,
    'Out_of_range': None}

HKcorr['6.5.X'] = HKcorr['6.3.0'].copy()

HKcorr['7.2.X'] = HKcorr['6.5.X'].copy()
ignore = list(map(HKcorr['7.2.X'].pop, ['fpgaSpwErr']))
HKcorr['7.2.X'].update({'dacExeErr': None,
                        'spwTmRspErrFlg': None,
                        'spwTcPktErrFlg': None,
                        'adcBsyErrFlg': None,
                        'CDPUClkLost': None})

HKcorr['7.5.X'] = HKcorr['7.2.X'].copy()
ignore = list(map(HKcorr['7.2.X'].pop, ['CDPUClkLost', 'CDPUClkSt']))
HKcorr['7.5.X'].update(dict(CDPUClkSw=None, pixelToErr=None, spwStatus=None, spwTmPtypeFlag=None,
                            TwoCmdsErr=None))

allstats = ['mean', 'std', 'min', 'max']


def loadHK_preQM(filename, elvis='5.7.07'):
    """Loads a HK file

    It only assumes a structure given by a HK keyword followed by a number of
    of tab-separated values (number not specified).
    Note that the length of the values arrays is variable (depends on length
    of exposure and HK sampling rate).

    :param filename: path to the file to be loaded, including the file itself

    :return: dictionary with pairs parameter:[values]

    """

    with open(filename) as f:
        lines = f.readlines()
        f.close()

    data = {}

    for line in lines:
        items = line.split()
        key = items[0]
        values = np.array(items[1:]).astype('float32')

        data[key] = values

    return data


def loadHK_QFMsingle(filename, elvis=context.elvis, validate=False, safe=False):
    """Loads a HK file

    Structure: tab separated columns, one per Keyword. First column is a
    timestamp, and there may be a variable number of rows (readings).

    :param filename: path to the file to be loaded, including the file itself
    :param elvis: "ELVIS" version

    :return: astropy table with pairs parameter:[values]

    """

    if safe:
        def wrapped_safe_open(path):
            return utils.safe_open_file(path, prefix='HK_')
        opener = wrapped_safe_open
    else:
        opener = open

    with opener(filename) as f:
        lines = f.readlines()
        f.close()
        table = ascii.read(lines)

    if validate:
        expectedkeys = allHK_keys[elvis]
        assert list(table.keys()) == expectedkeys, \
            'HK keys in %s not mathing expectation for ELVIS=%s' % \
                         (filename, elvis)

    return table


def mergeHK(HKList):
    """ """
    return vstack(HKList)

    #HK = copy.deepcopy(HKList[0])
    #nHK = len(HKList)

    # for ixHK in range(1,nHK):
    #    iHKdata = HKList[ixHK]

    # Row-by-Row appending of the "iHKdata" added catalog... is there a faster way?

    #    for iL in range(len(iHKdata)):
    #        HK.add_row(iHKdata[iL].as_void().tolist())

    # return HK


def loadHK_QFM(filename, elvis=context.elvis, validate=False, safe=False):
    """Loads a HK file, or list of HK files.

    Structure: astropy table. First column is a timestamp, and there may be
    a variable number of rows (readings).

    :param filename: path to the file to be loaded, including the file itself, or list of paths to HK files.
    :param elvis: "ELVIS" version

    :return: astropy table with pairs parameter:[values]

    """

    if isinstance(filename, str):
        return loadHK_QFMsingle(filename, elvis=elvis, validate=validate, safe=safe)
    elif isinstance(filename, list):
        HKs = [loadHK_QFMsingle(item, elvis=elvis, validate=validate, safe=safe)
               for item in filename]
        return mergeHK(HKs)


def iniHK_QFM(elvis=context.elvis, length=0):
    """ """
    columns = allHK_keys[elvis]

    dtypes = ['S8']  # TimeStamp
    for column in columns[1:]:
        dtypes.append('f')

    if length > 0:

        ColList = []
        for i, column in enumerate(columns):
            ColList.append(Column(name=column, dtype=dtypes[i], length=length))
        HK = Table(ColList)

    else:
        HK = Table(names=columns, dtype=dtypes)

    return HK


def filtervalues(values, key):
    """ """

    if ('HK_temp_top_CCD' in key) or ('HK_temp_bottom_CCD' in key):
        defvalue = -246.86
    elif ('HK_Temp1_RPSU' in key) or ('HK_Temp2_RPSU' in key):
        defvalue = -273.27
    elif (key == 'Video_TOP') or (key == 'Video_BOT'):
        defvalue = -273.15
    else:
        defvalue = 0.

    try:
        return values[np.where(values != defvalue)]
    except BaseException:
        return np.zeros_like(values) + np.nan


def synthHK(HK):
    """
    Synthetizes the values for each parameter in a HK dictionary into
    [mean,std,min,max].


    :param HK: a dictionary as those output by loadHK.

    :return: dictionary with pairs parameter:[mean,std,min,max]

    """

    synthHKdict = {}

    for key in list(HK.keys()):
        values = HK[key].data.copy()
        try:
            mean = values.mean()
        except BaseException:
            mean = np.nan
        try:
            std = values.std()
        except BaseException:
            std = np.nan
        try:
            minv = values.min()
        except BaseException:
            minv = np.nan
        try:
            maxv = values.max()
        except BaseException:
            maxv = np.nan

        synthHKdict[key] = np.array([mean, std, minv, maxv])

    return synthHKdict


def reportHK(HKs, key, reqstat='all'):
    """Returns (mean, std, min, max) for each keyword in a list of
    HK dictionaries (output from loadHK).

    :param HK: dictionary with HK data.
    :param key: HK key.
    :reqstat: what statistic to retrieve.

    """

    if reqstat != 'all':
        ixstat = allstats.index(reqstat)
    else:
        ixstat = np.where(allstats)

    vals = []

    for HK in HKs:
        vals.append(synthHK[HK][ixstat])

    vals = np.array(vals)

    return vals


def parseDTstr(DTstr):
    """ """

    date = DTstr[0:DTstr.index('D')]

    y2d = int(date[4:6])
    if y2d < 20:
        century = 2000
    else:
        century = 1900
    dd, MM, yy = int(date[0:2]), int(date[2:4]), y2d + century

    time = DTstr[DTstr.index('D') + 1:-1]
    hh, mm, ss = int(time[0:2]), int(time[2:4]), int(time[4:6])

    DTobject = datetime.datetime(yy, MM, dd, hh, mm, ss)

    return DTobject


def parseHKfname(HKfname):
    """Parses name of a HK file to retrieve OBSID, date and time, and ROE number.

    :parameter HKfname: name of HK file.
    :return: obsid,dtobj=datetime.datetime(yy,MM,dd,hh,mm,ss),ROE


    """

    stripHKfname = os.path.basename(HKfname)
    stripHKfname = os.path.splitext(stripHKfname)[0]

    items = stripHKfname.split( '_')

    assert items[0] == 'HK'  # just a basic check

    obsid = int(items[1])
    DT = items[2]
    ROE = items[3]

    DTobject = parseDTstr(DT)

    return obsid, DTobject, ROE


def parseHKfiles(HKlist, elvis=context.elvis):
    """

    :param HKlist: list of HK files (path+name).
    :param elvis: "ELVIS" version.

    :return: [obsids],[dtobjs],[tdeltasec],[HK_keys], [data(nfiles,nstats,nHKparams)]

    """

    nfiles = len(HKlist)
    HK_keys = allHK_keys[elvis]

    nkeys = len(HK_keys)

    ver = int(elvis.split('.')[0])
    subver = int(elvis.split('.')[1])

    if ver == 5 and subver <= 7:
        HKloader = loadHK_preQM
    else:
        HKloader = loadHK_QFM

    obsids = []
    dtobjs = []

    data = np.zeros((nfiles, 4, nkeys), dtype='float32')

    for ix, HKfname in enumerate(HKlist):

        obsid, dtobj, ROE = parseHKfname(HKfname)
        obsids.append(obsid)
        dtobjs.append(dtobj)

        HK = HKloader(HKfname, elvis=elvis)

        synthHKdict = synthHK(HK)

        for ik, key in enumerate(HK_keys):
            try:
                data[ix, :, ik] = synthHKdict[key]
            except BaseException:
                data[ix, :, ik] = np.nan  # HACK

    # reordering

    tsort = np.argsort(dtobjs)
    obsids = np.array(obsids)[tsort]
    dtobjs = np.array(dtobjs)[tsort]
    data = np.array(data)[tsort, :, :]

    tdeltasec = np.array([(item - dtobjs[0]).seconds for item in dtobjs])

    return obsids, dtobjs, tdeltasec, HK_keys, data


def format_date(x, pos=None):
    try:
        return pylab.num2date(x).strftime('%d/%b-%H:%M')
    except BaseException:
        pass
    # except: stop()

# ==============================================================================
# def HKplot(allHKdata,keylist,key,dtobjs,filename='',stat='mean'):
#     """Plots the values of a HK parameter as a function of time.
#
#     :param allHKdata: HKdata = [(nfiles,nstats,nHKparams)]
#     :param keylist: list with all HK keys.
#     :param key: selected key.
#     :param dtobjs: datetime objects time axis.
#     :param filename: file-name to store plot [empty string not to save].
#     :param stat: statistics to plot.
#
#     :return: None!!
#
#     """
#
#     ixkey = keylist.index(key)
#
#     ixstat = allstats.index(stat)
#
#     HKvals = allHKdata[:,ixstat,ixkey]
#
#     ylabel = 'V'
#
#     curr_keys = ['hk_vrclk_lo_roe','hk_i3.3v_dig_rpsu','hk_i1.5v_dig_rpsu',
#                  'hk_i28v_rpsu','hk_i+van_rpsu','hk_i-van_rpsu',
#                  'hk_ivclk_rpsu','hk_ivccd_rpsu','hk_vclk_roe']
#
#     if key.lower() in curr_keys:
#         ylabel = 'I[A]'
#
#     temp_keys = ['hk_temp_top_ccd1','hk_temp_bottom_ccd1',
#                  'hk_temp_top_ccd2','hk_temp_bottom_ccd2',
#                  'hk_temp_top_ccd3','hk_temp_bottom_ccd3',
#                  'hk_temp1_rpsu','hk_temp2_rpsu',
#                  'video_top','video_bot']
#
#     if key.lower() in temp_keys:
#         ylabel = 'deg [C]'
#
#     fig = plt.figure(figsize=(7,6))
#     ax1 = fig.add_subplot(111)
#     if not np.all(np.isnan(HKvals)):
#         ax1.plot(dtobjs,HKvals,'b-')
#     #ax1.plot(HKvals,'b')
#     #ax1.set_xlabel(r'$Time$')
#     ax1.set_ylabel(r'$%s$' % ylabel)
#     ntitle = key.replace('_','\_')
#     ax1.set_title(ntitle)
#
#     for item in [ax1.title, ax1.xaxis.label, ax1.yaxis.label]:
#         item.set_fontsize(18)
#
#     ylim = ax1.get_ylim()
#     dy = ylim[1]-ylim[0]
#     if dy < 0.1:
#         ymid = (ylim[1]+ylim[0])/2.
#         ylim = (ymid-0.05,ymid+0.05)
#
#     ax1.set_ylim(ylim)
#
#     #if key == 'HK_temp_top_CCD1': stop()
#
#     #if dtobjs[0] != 0:
#     plt.gca().xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(format_date))
#     fig.autofmt_xdate()
#
#
#     plt.tight_layout(rect=[0, 0, 1, 1])
#
#     if filename== '':
#         plt.show()
#     else:
#         plt.savefig(filename)
#
#
#     plt.close()
#
# ==============================================================================


def _ax_render_HK(ax, x, y, HKlims, HKkey, fontsize=10):
    """ """

    max_xticks = 6
    xloc = plt.MaxNLocator(max_xticks)

    ax.clear()

    if np.any(np.isnan(y)):
        yp = y.copy()
        yp[np.isnan(y)] = 0
        ax.plot(x, yp)
        ax.plot(x[np.isnan(y)], yp[np.isnan(y)], 'ro')
    else:
        ax.plot(x, y, 'b.', ms=3)

    if len(HKlims) == 1:
        ax.axhline(y=HKlims[0], ls='--', lw=2, color='r')
    elif len(HKlims) == 2:
        for ik in range(2):
            ax.axhline(y=HKlims[ik], ls='--', lw=2, color='r')

    HKtitle = '$%s$' % HKkey.replace('_', '\_')
    ax.set_title(HKtitle)

    try:
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                     ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(fontsize)
    except BaseException:
        pass

    ax.xaxis.set_major_locator(xloc)

    return ax


def doHKSinglePlot(dtobjs, HK, HKkey, ylabel='V', HKlims=[], filename='', fontsize=10):
    """Plots the values of a HK parameter as a function of time.

    :param dtobjs: datetime objects time axis.
    :param HK: HK values (array)
    :param HKkey:
    :param ylabel:
    :param HKlims:
    :param filename: file-name to store plot [empty string not to save].

    :return: None!!

    """

    fig = plt.figure(figsize=(7, 6))
    ax1 = fig.add_subplot(111)

    ax1 = _ax_render_HK(ax1, dtobjs, HK, HKlims, HKkey, fontsize=fontsize)

    plt.gca().xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(format_date))
    fig.autofmt_xdate()

    plt.tight_layout(rect=[0, 0, 1, 1])

    if filename == '':
        plt.show()
    else:
        plt.savefig(filename)

    plt.close()


def check_HK_vs_command(HKKeys, dd, limits='P', elvis=context.elvis):
    """
    Returns report on HK parameters, in DataDict (dd), comparing inputs (commanded)
    vs. output (HK data).

    HK Keys which do not correspond to commanded voltages always return 'True'.

    :param HKKeys: list of HK parameters, as named in HK files (without HK_ suffix)
    :param dd: DataDict object
    :param limits: type of limits to use, either "P" (Performance) or "S" (Safe)
    :param elvis: ELVIS version to find correspondence between HK key and
                 Exposure Log input (commanded voltage).
    :returns report: dictionary with pairs of HK-key : Bool.
                 True = All values are within limits, referred to commanded value.
                 False = At least one value is outside limits, referred to commanded value.

    """

    #raise NotImplementedError

    report = dict()

    for HKKey in HKKeys:

        HKlim = HKlims[elvis][limits][HKKey]
        limtype = HKlim[0]
        limval = HKlim[1:]

        ELKey = HKcorr[elvis][HKKey]
        if ELKey is None:
            report[HKKey] = True
            continue

        ELdata = dd.mx[ELKey][:]
        HKdata = dd.mx['HK_%s' % HKKey][:, np.newaxis]

        if limtype == 'R':
            test = HKdata - ELdata
            testBool = ((test <= limval[0]) | (test >= limval[1]))
            report[HKKey] = not np.any(testBool)
        elif limtype == 'A':
            testBool = ((HKdata <= limval[0]) | (test >= limval[1]))
            report[HKKey] = not np.any(testBool)
        elif limtype == 'I':
            testBool = HKdata != limval[0]
            report[HKKey] = not np.all(test)

        # print HKKey, limtype, test.mean(), report[HKKey] # TESTS

    return report


def check_HK_abs(HKKeys, dd, limits='S', elvis=context.elvis):
    """

    Returns report on HK parameters, in DataDict (dd), compared to absolute
    limits.

    HK Keys which have "relative" limits, always return False.

    :param HKKeys: list of HK parameters, as named in HK files (without HK_ suffix)
    :param dd: DataDict object
    :param limits: type of limits to use, either "P" (Performance) or "S" (Safe)
    :param elvis: ELVIS version to find correspondence between HK key and
                 Exposure Log input (commanded voltage).
    :returns report: dictionary with pairs of HK-key : Bool.
                 True = All values for given key are within limits.
                 False = At least one value for given key is outside limits.


    """

    report = dict()

    for HKKey in HKKeys:

        HKlim = HKlims[elvis][limits][HKKey]
        limtype = HKlim[0]
        limval = HKlim[1:]

        HKdata = dd.mx['HK_%s' % HKKey][:, np.newaxis]

        if limtype == 'R':
            report[HKKey] = False
        elif limtype == 'A':
            testBool = (HKdata <= limval[0]) | (HKdata >= limval[1])
            report[HKKey] = not np.any(testBool)
        elif limtype == 'I':
            testBool = HKdata != limval[0]
            report[HKKey] = not np.all(test)

        # print HKKey, limtype, report[HKKey] # TESTS

    return report


def test():

    datapath = 'data/08_Mar/'
    HKfile = os.path.join(datapath, 'HK_281_080316D174300T_ROE1.txt')

    #HK = loadHK(HKfile)
    #HKbrief = synthHK(HK)

    HKlist = glob('%s/HK_*ROE1.txt' % datapath)

    obsids, dtobjs, tdeltasec, HK_keys, allHKdata = parseHKfiles(HKlist)
    # [allHKdata] = [n_obsids, n_stats, n_hk]
    tdeltahour = tdeltasec / 3600.
    HK_temp_top_CCD2 = allHKdata[:, 0, HK_keys.index('HK_temp_top_CCD2')]
    HK_temp_bot_CCD2 = allHKdata[:, 0, HK_keys.index('HK_temp_bottom_CCD2')]

    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax1.plot(tdeltahour, HK_temp_bot_CCD2, 'b-', label='bottom')
    ax1.plot(tdeltahour, HK_temp_top_CCD2, 'r-', label='top')
    ax1.set_xlabel(r'$\Delta T[hour]$')
    ax1.set_ylabel('Temps CCD2 [C]')
    ax1.set_ylim([-120., -40.])

    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles, labels)

    ax2 = fig.add_subplot(122)
    ax2.plot(tdeltahour, HK_temp_bot_CCD2, 'b-')
    ax2.plot(tdeltahour, HK_temp_top_CCD2, 'r-')
    ax2.set_xlabel(r'$\Delta T[hour]$')
    ax2.set_ylabel('Temps CCD2 [C]')
    # ax2.set_ylim([-105,-95])

    fig.suptitle('08.March.2016 - Functionals')
    plt.show()

    for key in HK_keys:
        HKplot(allHKdata, HK_keys, key, tdeltahour)


if __name__ == '__main__':

    test()
