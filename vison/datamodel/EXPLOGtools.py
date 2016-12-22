#! /home/raf/SOFTWARE/anaconda2/envs/VISSIM/bin/python

# IMPORT STUFF
from astropy.io import ascii
import os
from copy import copy

from pdb import set_trace as stop

# END IMPORT

columnlist = dict(v5702=['ObsID','File_name','CCD','ROE','DATE','PROGRAM','TEST','CCD_SN','BUNIT',
    'Operator','Lab_ver','Con_file','Exptime','N_P_high','Chrg_inj','On_cycle',
    'Off_cycl','Rpeat_cy','pls_leng','pls_del','Trappump','TP_Ser_S',
    'TP_Ver_S','TP_DW_V','TP_DW_H','TOI_flu','TOI _pump','TOI_read','Invflshp',
    'Invflush','Flushes','Vstart','Vend','Lim_scan','Ovrscn_H','Ovrscn_V',
    'V_clocks','S_clocks','CLK_ROE1','CLK_ROE2','ROE1_On','ROE2_On','FPGA_ver','EGSE_ver',
    'M_Steps','M_St_Sze','Wavelength','Lght_int','Chmb_pre','CCD_SNX','ROE_SN',
    'CalScrpt','R1CCD1TT','R1CCD1TB','R1CCD2TT','R1CCD2TB','R1CCD3TT',
    'R1CCD3TB','R2CCD1TT','R2CCD1TB','R2CCD2TT','R2CCD2TB',
    'R2CCD3TT','R2CCD3TB','IDL_V','IDH_V','IG1_V','IG2_V','ODCCD1_V',
    'RDCCD1_V','ODCCD2_V','RDCCD2_V','ODCCD3_V','RDCCD3_V'],
    v5704=['ObsID','File_name','CCD','ROE','DATE','PROGRAM','TEST','CCD1_SN',
              'CCD2_SN','CCD3_SN','BUNIT','Operator','Lab_ver','Con_file',
              'Exptime','Flsh-Rdout_e_time','C.Inj-Rdout_e_time','N_P_high',
              'Chrg_inj','On_cycle','Off_cycl','Rpeat_cy','pls_leng','pls_del',
              'Trappump','TP_Ser_S','TP_Ver_S','TP_DW_V','TP_DW_H','TOI_flu',
              'TOI _pump','TOI_read','Invflshp','Invflush','Flushes','Vstart',
              'Vend','	Lim_scan','Ovrscn_H','Ovrscn_V','V_clocks','S_clocks',
              'CLK_ROE1','CLK_ROE2','ROE1_On','ROE2_On','FPGA_ver','EGSE_ver',
              'M_Steps	','M_St_Sze','Wavelength','Lght_int','Chmb_pre','ROE_SN',
              'CalScrpt','R1CCD1TT','R1CCD1TB','R1CCD2TT','R1CCD2TB','R1CCD3TT',
              'R1CCD3TB','R2CCD1TT','R2CCD1TB','R2CCD2TT','R2CCD2TB',
              'R2CCD3TT','R2CCD3TB','IDL_V','IDH_V','IG1_V','IG2_V','ODCCD1_V',
              'RDCCD1_V','ODCCD2_V','RDCCD2_V','ODCCD3_V','RDCCD3_V'])

def loadExpLog(expfile,elvis='v5704'):
    """ """
    
    explog = ascii.read(expfile,data_start=1,delimiter='\t',guess=False,\
                        names=columnlist[elvis],format='no_header')
    return explog
    
def mergeExpLogs(explogList):
    """ """
    colnames = explogList[0].colnames
    
    nexplogs=len(explogList)
    
    assert nexplogs >=2
    
    explog = copy(explogList[0])
    
    for iexp in range(1,nexplogs):
        iexplog = explogList[iexp]
        
        for iL in range(len(iexplog)):
            explog.add_row(iexplog[iL].as_void().tolist())
    
    
    return explog


def test():
    """ """
    expfile1 = os.path.join('14_Mar_16L','EXP_LOG_140316.txt')
    expfile2 = os.path.join('15_Mar_16L','EXP_LOG_150316.txt')
    
    explog1 = loadExpLog(expfile1)
    explog2 = loadExpLog(expfile2)
    
    explog = mergeExpLogs([explog1,explog2])
    
    stop()
    
if __name__ == '__main__':
    """Tests"""
    
    test()
    
