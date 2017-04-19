#! /home/raf/SOFTWARE/anaconda2/envs/VISSIM/bin/python

# IMPORT STUFF
from pdb import set_trace as stop


from astropy.io import ascii
from astropy.table import Table, Column
import os
from copy import copy
import numpy as np


# END IMPORT

columnlist = {'5.7.02':['ObsID','File_name','CCD','ROE','DATE','PROGRAM','TEST','CCD_SN','BUNIT',
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
    '5.7.04':['ObsID','File_name','CCD','ROE','DATE','PROGRAM','TEST','CCD1_SN',
              'CCD2_SN','CCD3_SN','BUNIT','Operator','Lab_ver','Con_file',
              'Exptime','Flsh-Rdout_e_time','C.Inj-Rdout_e_time','N_P_high',
              'Chrg_inj','On_cycle','Off_cycl','Rpeat_cy','pls_leng','pls_del',
              'Trappump','TP_Ser_S','TP_Ver_S','TP_DW_V','TP_DW_H','TOI_flu',
              'TOI_pump','TOI_read','Invflshp','Invflush','Flushes','Vstart',
              'Vend','Lim_scan','Ovrscn_H','Ovrscn_V','V_clocks','S_clocks',
              'CLK_ROE1','CLK_ROE2','ROE1_On','ROE2_On','FPGA_ver','EGSE_ver',
              'M_Steps','M_St_Sze','Wavelength','Lght_int','Chmb_pre','ROE_SN',
              'CalScrpt','R1CCD1TT','R1CCD1TB','R1CCD2TT','R1CCD2TB','R1CCD3TT',
              'R1CCD3TB','R2CCD1TT','R2CCD1TB','R2CCD2TT','R2CCD2TB',
              'R2CCD3TT','R2CCD3TB','IDL_V','IDH_V','IG1_V','IG2_V','ODCCD1_V',
              'RDCCD1_V','ODCCD2_V','RDCCD2_V','ODCCD3_V','RDCCD3_V']}

columnlist['5.7.09'] = columnlist['5.7.04']

columnlist['6.0.0'] = ['ObsID','File_name','CCD','ROE','DATE','PROGRAM',
'TEST','CCD1_SN','CCD2_SN','CCD3_SN','BUNIT','Operator','Lab_ver','Con_file',
'Exptime','Flsh-Rdout_e_time','C.Inj-Rdout_e_time','N_P_high','Chrg_inj',
'On_cycle','Off_cycl','Rpeat_cy','pls_leng','pls_del','SerRdDel','Trappump',
'TP_Ser_S','TP_Ver_S','TP_DW_V','TP_DW_H','TOI_flsh','TOI_pump','TOI_read',
'TOI_CInj','Invflshp','Invflush','Flushes','Vstart','Vend','Ovrscn_H',
'CLK_ROE','CnvStart','SumWell','IniSweep','SPW_clk','FPGA_ver','EGSE_ver',
'M_Steps','M_St_Sze','Wavelength','Mirr_pos','Chmb_pre','RPSU_SN','ROE_SN','CalScrpt',
'R1CCD1TT','R1CCD1TB','R1CCD2TT','R1CCD2TB','R1CCD3TT','R1CCD3TB','IDL_V',
'IDH_V','IG1_T1_V','IG1_T2_V','IG1_T3_V','IG1_B1_V','IG1_B2_V','IG1_B3_V',
'IG2_T_V','IG2_B_V','OD_T1_V','OD_T2_V','OD_T3_V','OD_B1_V','OD_B2_V',
'OD_B3_V','RD_T_V','RD_B_V']


dtypes_dict = {'ObsID':'i','File_name':'S40','CCD':'S4','ROE':'S5','DATE':'S20',
              'PROGRAM':'S12','TEST':'S12','CCD1_SN':'S20','CCD2_SN':'S20',
              'CCD3_SN':'S20','BUNIT':'S3','Operator':'S3','Lab_ver':'S8',
              'Con_file':'S30','Exptime':'f',
              'Flsh-Rdout_e_time':'f','C.Inj-Rdout_e_time':'f',
              'N_P_high':'S6','Chrg_inj':'i', 'On_cycle':'i','Off_cycl':'i',
              'Rpeat_cy':'i','pls_leng':'f','pls_del':'f','SerRdDel':'i','Trappump':'S10', 
              'TP_Ser_S':'i','TP_Ver_S':'i','TP_DW_V':'i','TP_DW_H':'i','TOI_flsh':'i','TOI_flu':'i',
              'TOI_pump':'i','TOI_read':'i','TOI_CInj':'f','Invflshp':'i','Invflush':'i',
              'Flushes':'i','Vstart':'i','Vend':'i','Ovrscn_H':'i','CLK_ROE':'S10',
              'CnvStart':'i','SumWell':'i','IniSweep':'i','SPW_clk':'i','FPGA_ver':'S10','EGSE_ver':'S10',
              'M_Steps':'f','M_St_Sze':'f','Wavelength':'i','Mirr_pos':'f','RPSU_SN':'S10','ROE_SN':'S10',
              'CalScrpt':'S30','R1CCD1TT':'f','R1CCD1TB':'f','R1CCD2TT':'f','R1CCD2TB':'f','R1CCD3TT':'f',
              'R1CCD3TB':'f','IDL_V':'f','IDH_V':'f','IG1_T1_V':'f','IG1_T2_V':'f','IG1_T3_V':'f',
              'IG1_B1_V':'f','IG1_B2_V':'f','IG1_B3_V':'f','IG2_T_V':'f','IG2_B_V':'f',
              'OD_T1_V':'f','OD_T2_V':'f','OD_T3_V':'f','OD_B1_V':'f','OD_B2_V':'f',
              'OD_B3_V':'f','RD_T_V':'f','RD_B_V':'f','TOI_flsh':'i',
              'Lim_scan':'i','Ovrscn_V':'i','V_clocks':'S10','S_clocks':'S10',
              'CLK_ROE1':'S10','CLK_ROE2':'S10','ROE1_On':'i','ROE2_On':'i',
              'Lght_int':'f','Chmb_pre':'f','R2CCD1TT':'f','R2CCD1TB':'f','R2CCD2TT':'f',
              'R2CCD2TB':'f','R2CCD3TT':'f','R2CCD3TB':'f',
              'IG1_V':'f','IG2_V':'f','ODCCD1_V':'f',
              'RDCCD1_V':'f','ODCCD2_V':'f','RDCCD2_V':'f','ODCCD3_V':'f','RDCCD3_V':'f'}


def iniExplog(elvis):
    """ """
    columns = columnlist[elvis]    
    
    dtypes = []
    for column in columns:
        dtypes.append(dtypes_dict[column])
    
    explog = Table(names=columns,dtype=dtypes)
    
    return explog
    


def loadExpLog(expfile,elvis='5.7.04'):
    """Loads an Exposure Log from file."""
    explog = ascii.read(expfile,data_start=1,delimiter='\t',guess=False,\
                        names=columnlist[elvis],format='no_header')
    return explog
    

def mergeExpLogs(explogList,addpedigree=False):
    """Merges explog objects in a list."""
    
    nexplogs=len(explogList)
    
    #assert nexplogs >=2
    
    explog = copy(explogList[0])
    colnames = explog.colnames
    
    if addpedigree:
        explog['explognumber'] = np.zeros(len(explog),dtype='int32')
    
    for iexp in range(1,nexplogs):
        iexplog = explogList[iexp]
        
        # Check format compatibility of columns
        
        for ic,colname in enumerate(colnames):
            
            if iexplog[colname].dtype.kind == 'S' and not explog[colname].dtype.kind == 'S':
                explog[colname] = explog[colname].astype(str)
        
        if addpedigree:
            iexplog['explognumber'] = np.ones(len(iexplog),dtype='int32') * iexp
        
        # Row-by-Row appending of the "iexp" added catalog
        
        for iL in range(len(iexplog)):
            explog.add_row(iexplog[iL].as_void().tolist())
    
    
    return explog

class ExpLogClass():
    """ """
    
    def __init__(self,elvis='5.7.04'):
        """ """
        
        self.elvis = elvis
        
    def loadfromfile(self,expfile,elvis='5.7.04'):
        
        self.explog = loadExpLog(expfile,self.elvis)
    
    def writeto(self,outfile):
        """ """
        
        self.explog.write(outfile,format='ascii')
    
    def iniExplog(self):
        """ """
        
        self.explog = iniExplog(self.elvis)
    
        return None
    
    
    def addRow(self):
        """ """
    


def test():
    """This Tests needs UPDATE (for data access and probably data format)"""
    expfile1 = os.path.join('14_Mar_16L','EXP_LOG_140316.txt')
    expfile2 = os.path.join('15_Mar_16L','EXP_LOG_150316.txt')
    
    explog1 = loadExpLog(expfile1)
    explog2 = loadExpLog(expfile2)
    
    explog = mergeExpLogs([explog1,explog2])
    
    stop()
    
if __name__ == '__main__':
    """Tests"""
    
    test()
    
