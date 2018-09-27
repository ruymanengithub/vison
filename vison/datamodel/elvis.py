#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

ELVIS variables dictionaries.

Created on Fri Sep 22 12:04:09 2017

:author: Ruyman Azzollini

"""

script_keys = ['frames',
               'program', 'test', 'IDL', 'IDH', 'IG1_1_T', 'IG1_2_T', 'IG1_3_T', 'IG1_1_B',
               'IG1_2_B', 'IG1_3_B', 'IG2_T', 'IG2_B', 'OD_1_T', 'OD_1_B', 'OD_2_T', 'OD_2_B',
               'OD_3_T', 'OD_3_B', 'RD_T', 'RD_B', 'IPHI1', 'IPHI2', 'IPHI3', 'IPHI4', 'rdmode',
               'flushes', 'siflsh', 'siflsh_p', 'swellw', 'swelldly', 'inisweep', 'vstart', 'vend',
               'toi_fl', 'toi_tp', 'toi_ro', 'toi_ch', 'chinj', 'chinj_on', 'chinj_of', 'id_wid',
               'id_dly', 'chin_dly', 's_tpump', 's_tpmod', 'v_tpump', 'v_tpmod', 's_tp_cnt', 'v_tp_cnt',
               'dwell_v', 'dwell_s', 'exptime', 'shuttr', 'e_shuttr', 'pos_mirr', 'wave', 'motr',
               'motr_cnt', 'motr_siz', 'source', 'operator', 'sn_ccd1', 'sn_ccd2', 'sn_ccd3', 'sn_roe',
               'sn_rpsu', 'comments']


explog_keys = ['ObsID',
               'File_name', 'CCD', 'ROE', 'date', 'program', 'test', 'calscrpt', 'sn_ccd1',
               'sn_ccd2', 'sn_ccd3', 'sn_roe', 'sn_rpsu', 'BUNIT', 'operator', 'con_file',
               'exptime', 'fl_rdout', 'ci_rdout', 'IPHI', 'vstart', 'vend', 'rdmode', 'flushes',
               'siflsh', 'siflsh_p', 'swellw', 'swelldly', 'inisweep', 'spw_clk', 'chinj', 'chinj_on',
               'chinj_of', 'id_wid', 'id_dly', 'chin_dly', 'v_tpump', 's_tpump', 'v_tp_mod', 's_tp_mod',
               'v_tp_cnt', 's_tp_cnt', 'dwell_v', 'dwell_s', 'toi_fl', 'toi_tp', 'toi_ro', 'toi_ch',
               'fpga_ver', 'egse_ver', 'motr', 'motr_cnt', 'motr_siz', 'source', 'wave', 'pos_mirr',
               'R1C1_TT', 'R1C1_TB', 'R1C2_TT', 'R1C2_TB', 'R1C3_TT', 'R1C3_TB', 'IDL', 'IDH',
               'IG1_1_T', 'IG1_2_T', 'IG1_3_T', 'IG1_1_B', 'IG1_2_B', 'IG1_3_B', 'IG2_T',
               'IG2_B', 'OD_1_T', 'OD_2_T', 'OD_3_T', 'OD_1_B',
               'OD_2_B', 'OD_3_B', 'RD_T', 'RD_B']


FITS_keys = ['XTENSION', 'BITPIX', 'NAXIS', 'NAXIS1', 'NAXIS2', 'PCOUNT', 'GCOUNT',
             'BZERO', 'BSCALE', 'EXTNAME', 'BUNIT', 'PROGRAM', 'OBJECT', 'OBSID', 'OPERATOR',
             'FULLPATH', 'LAB_VER', 'CON_FILE', 'DATE', 'EXPTIME', 'FL_RDOUT', 'CI_RDOUT', 'IPHI',
             'CHINJ', 'CHINJ_ON', 'CHINJ_OF', 'ID_WID', 'ID_DLY', 'chin_dly', 'V_TPUMP',
             'S_TPUMP', 'V_TP_MOD', 'S_TP_MOD', 'V_TP_CNT', 'S_TP_CNT', 'DWELL_V', 'DWELL_S',
             'TOI_FL', 'TOI_TP', 'TOI_RO', 'TOI_CH', 'SIFLSH', 'SIFLSH_P', 'FLUSHES', 'VSTART',
             'VEND', 'RDMODE', 'SWELLW', 'SWELLDLY', 'INISWEEP', 'SPW_CLK', 'FPGA_VER', 'EGSE_VER',
             'MOTR', 'MOTR_CNT', 'MOTR_SIZ', 'SOURCE', 'WAVE', 'POS_MIRR', 'SN_CCD1', 'SN_CCD2',
             'SN_CCD3', 'SN_ROE', 'SN_RPSU', 'CALSCRPT', 'COMMENTS', 'R1C1_TT', 'R1C1_TB', 'R1C2_TT',
             'R1C2_TB', 'R1C3_TT', 'R1C3_TT', 'IDL', 'IDH', 'IG1_1_T', 'IG1_2_T', 'IG1_3_T',
             'IG1_1_B', 'IG1_2_B', 'IG1_3_B', 'IG2_T', 'IG2_B', 'OD_1_T', 'OD_2_T', 'OD_3_T',
             'OD_1_B', 'OD_2_B', 'OD_3_B', 'RD_T', 'RD_B']
