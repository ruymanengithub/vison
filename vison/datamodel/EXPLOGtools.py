#! /home/raf/SOFTWARE/anaconda2/envs/VISSIM/bin/python

# IMPORT STUFF
from pdb import set_trace as stop


from astropy.io import ascii
from astropy.table import Table, Column
import os
from copy import copy
import numpy as np
from vison.support import vistime
from vison.pipe import lib as pilib

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
'On_cycle','Off_cycl','pls_leng','pls_del','SerRdDel','Trappump',
'TP_Ser_S','TP_Ver_S','TP_DW_V','TP_DW_H','TOI_flsh','TOI_pump','TOI_read',
'TOI_CInj','Invflshp','Invflush','Flushes','Vstart','Vend','Ovrscn_H',
'CLK_ROE','CnvStart','SumWell','IniSweep','SPW_clk','FPGA_ver','EGSE_ver',
'M_Steps','M_St_Sze','Wavelength','Mirr_pos','Chmb_pre','RPSU_SN','ROE_SN','CalScrpt',
'R1CCD1TT','R1CCD1TB','R1CCD2TT','R1CCD2TB','R1CCD3TT','R1CCD3TB','IDL_V',
'IDH_V','IG1_T1_V','IG1_T2_V','IG1_T3_V','IG1_B1_V','IG1_B2_V','IG1_B3_V',
'IG2_T_V','IG2_B_V','OD_T1_V','OD_T2_V','OD_T3_V','OD_B1_V','OD_B2_V',
'OD_B3_V','RD_T_V','RD_B_V']
columnlist['6.1.0'] = columnlist['6.0.0']

columnlist['6.3.0'] = ['ObsID','File_name','CCD','ROE','date','program',
'test','calscrpt','sn_ccd1','sn_ccd2','sn_ccd3','sn_roe','sn_rpsu',
'BUNIT','operator','con_file','exptime','fl_rdout','ci_rdout','IPHI',
'vstart','vend','rdmode','flushes','siflsh','siflsh_p','swellw','swelldly',
'inisweep','spw_clk','chinj','chinj_on','chinj_of','id_wid','id_dly','chin_dly',
'v_tpump','s_tpump','v_tpmod','s_tpmod','v_tp_cnt','s_tp_cnt',
'dwell_v','dwell_s',
'toi_fl','toi_tp','toi_ro','toi_ch',
'fpga_ver','egse_ver',
'motr_on','motr_cnt','motr_siz','source','wave','mirr_on','mirr_pos',
'R1C1_TT','R1C1_TB','R1C2_TT','R1C2_TB','R1C3_TT','R1C3_TB','IDL','IDH',
'IG1_1_T','IG1_2_T','IG1_3_T','IG1_1_B','IG1_2_B','IG1_3_B','IG2_T','IG2_B',
'OD_1_T','OD_2_T','OD_3_T','OD_1_B','OD_2_B','OD_3_B','RD_T','RD_B']

columnlist['6.5.X'] = ['ObsID','File_name','CCD','ROE','date','program',
'test','calscrpt','sn_ccd1','sn_ccd2','sn_ccd3','sn_roe','sn_rpsu',
'BUNIT','operator','con_file','exptime','fl_rdout','ci_rdout','IPHI',
'vstart','vend','rdmode','flushes','siflsh','siflsh_p','swellw','swelldly',
'inisweep','cdpu_clk','chinj','chinj_on','chinj_of','id_wid','id_dly','chin_dly',
'v_tpump','s_tpump','v_tp_mod','s_tp_mod','v_tp_cnt','s_tp_cnt',
'dwell_v','dwell_s',
'toi_fl','toi_tp','toi_ro','toi_ch',
'fpga_ver','egse_ver',
'motr','motr_cnt','motr_siz','source','wave','mirr_on','mirr_pos',
'R1C1_TT','R1C1_TB','R1C2_TT','R1C2_TB','R1C3_TT','R1C3_TB','IDL','IDH',
'IG1_1_T','IG1_2_T','IG1_3_T','IG1_1_B','IG1_2_B','IG1_3_B','IG2_T','IG2_B',
'OD_1_T','OD_2_T','OD_3_T','OD_1_B','OD_2_B','OD_3_B','RD_T','RD_B']

dtypes_dict = {'ObsID':'i','File_name':'S40','CCD':'S4','ROE':'S5','DATE':'S20',
              'PROGRAM':'S12','TEST':'S12','CCD1_SN':'S20','CCD2_SN':'S20',
              'CCD3_SN':'S20','BUNIT':'S3','Operator':'S3','Lab_ver':'S8',
              'Con_file':'S40','Exptime':'f',
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
              'CLK_ROE':'S10','CLK_ROE1':'S10','CLK_ROE2':'S10','ROE1_On':'i','ROE2_On':'i',
              'Lght_int':'f','Chmb_pre':'f','R2CCD1TT':'f','R2CCD1TB':'f','R2CCD2TT':'f',
              'R2CCD2TB':'f','R2CCD3TT':'f','R2CCD3TB':'f',
              'IG1_V':'f','IG2_V':'f','ODCCD1_V':'f',
              'RDCCD1_V':'f','ODCCD2_V':'f','RDCCD2_V':'f','ODCCD3_V':'f','RDCCD3_V':'f',
              'date':'S20','program':'S12','test':'S12','calscrpt':'S30',
              'sn_ccd1':'S20','sn_ccd2':'S20','sn_ccd3':'S20',
              'sn_roe':'S20','sn_rpsu':'S20',
              'BUNIT':'S3','operator':'S3','con_file':'S40','exptime':'f',
              'fl_rdout':'f','ci_rdout':'f','IPHI':'S6',
              'vstart':'i','vend':'i','rdmode':'S12','flushes':'i','siflsh':'i',
              'siflsh_p':'f','swellw':'f','swelldly':'f',
              'inisweep':'i','spw_clk':'i','cdpu_clk':'i',
              'chinj':'i','chinj_on':'i',
              'chinj_of':'i','id_wid':'f','id_dly':'f','chin_dly':'i',
              'v_tpump':'i','s_tpump':'i','v_tpmod':'i','s_tpmod':'i',
              'v_tp_mod':'i','s_tp_mod':'i',
              'v_tp_cnt':'i','s_tp_cnt':'i',
              'dwell_v':'f','dwell_s':'f',
              'toi_fl':'i','toi_tp':'i','toi_ro':'i','toi_ch':'i',
              'fpga_ver':'S20','egse_ver':'S20',
              'motr_on':'i','motr':'i',
              'motr_cnt':'i','motr_siz':'f','source':'S6','wave':'i',
              'mirr_on':'i','mirr_pos':'f',
              'R1C1_TT':'f','R1C1_TB':'f','R1C2_TT':'f','R1C2_TB':'f',
              'R1C3_TT':'f','R1C3_TB':'f','IDL':'f','IDH':'f',
              'IG1_1_T':'f','IG1_2_T':'f','IG1_3_T':'f','IG1_1_B':'f',
              'IG1_2_B':'f','IG1_3_B':'f','IG2_T':'f','IG2_B':'f',
              'OD_1_T':'f','OD_2_T':'f','OD_3_T':'f','OD_1_B':'f',
              'OD_2_B':'f','OD_3_B':'f','RD_T':'f','RD_B':'f'}

#script_keys_cross = dict(
#              PROGRAM='program',TEST='test',CCD1_SN='sn_ccd1',CCD2_SN='sn_ccd2',
#              CCD3_SN='sn_ccd3',Operator='operator',Exptime='exptime',
#              Chrg_inj='chinj', On_cycle='chinj_rows_on',Off_cycl='chinj_rows_off',
#              pls_leng='id_width',pls_del='id_delay',SerRdDel='chinj_ser_wait',
#              Trappump='tpump', 
#              TP_Ser_S='ser_shuffles',TP_Ver_S='ver_shuffles',
#              TP_DW_V='dwell_v',TP_DW_H='dwell_h',TOI_flsh='toi_flush',TOI_flu='toi_flush',
#              TOI_pump='toi_tpump',TOI_read='toi_rdout',TOI_CInj='toi_chinj',
#              Invflshp='sinvflushp',Invflush='sinvflush',
#              Flushes='flushes',Vstart='vstart',Vend='vend',Ovrscn_H='add_h_overscan',
#              SumWell='sumwell',IniSweep='inisweep',
#              M_Steps='matrix_size',M_St_Sze='step_size',Wavelength='wavelength',
#              Mirr_pos='pos_cal_mirror',RPSU_SN='sn_rpsu',ROE_SN='sn_roe',
#              #IDL_V='f',IDH_V='f',IG1_T1_V='f',IG1_T2_V='f',IG1_T3_V='f',
#              #IG1_B1_V='f',IG1_B2_V='f',IG1_B3_V='f',IG2_T_V='f',IG2_B_V='f',
#              #OD_T1_V='f',OD_T2_V='f',OD_T3_V='f',OD_B1_V='f',OD_B2_V='f',
#              #OD_B3_V='f',RD_T_V='f',RD_B_V='f',TOI_flsh='i',
#              #Lim_scan='i',
#              Ovrscn_V='add_v_overscan',
#              #V_clocks='S10',S_clocks='S10',
#              CLK_ROE='readmode_1',CLK_ROE1='readmode_1',CLK_ROE2='readmode_2',
#              #ROE1_On='i',ROE2_On='i',
#              #Lght_int='f',Chmb_pre='f',R2CCD1TT='f',R2CCD1TB='f',R2CCD2TT='f',
#              #R2CCD2TB='f',R2CCD3TT='f',R2CCD3TB='f',
#              #IG1_V='f',IG2_V='f',ODCCD1_V='f',
#              #RDCCD1_V='f',ODCCD2_V='f',RDCCD2_V='f',ODCCD3_V='f',RDCCD3_V='f') 
#              )


def iniExplog(elvis):
    """ """
    columns = columnlist[elvis]    
    
    dtypes = []
    for column in columns:
        dtypes.append(dtypes_dict[column])
    
    explog = Table(names=columns,dtype=dtypes)
    
    return explog
    


def loadExpLog(expfile,elvis=pilib.elvis):
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
    
    def __init__(self,elvis=pilib.elvis):
        """ """
        self.elvis = elvis
        self.infile = ''
        
    def loadfromfile(self,expfile):
        
        self.infile = expfile
        self.explog = loadExpLog(expfile,self.elvis)
    
    def writeto(self,outfile):
        """ """
        self.explog.write(outfile,format='ascii')
    
    def iniExplog(self):
        """ """
        self.explog = iniExplog(self.elvis)
    
    
    def addRow(self,row):
        """ """
        self.explog.add_row(row.as_void().tolist())
    
    def summary(self):
        """ """
        try: nentries = len(self.explog)
        except AttributeError: 
            return None
        ObsID = self.explog['ObsID'].copy()
        nObsID = len(np.unique(ObsID))
        ObsIDextrema = np.min(ObsID),np.max(ObsID)
        
        tests = []
        for test in self.explog['test']: 
            if test not in tests: 
                tests.append(test)

        ntests = len(tests)
        
        durations = []
        for test in tests:
            ixsel = np.where(self.explog['test'] == test)
            ixlast = ixsel[0][-1]
            ixfirst = max(0,ixsel[0][0]-1)
            subdates = (self.explog['date'][ixfirst],self.explog['date'][ixlast])
            subdts = map(vistime.get_dtobj,subdates)
            dtmin = (subdts[-1]-subdts[0]).seconds / 60.
            durations.append(dtmin)
        
        totdurh = np.sum(durations) / 60.
        
        print '\nSummary: %s' % self.infile
        print 'Nr. of entries = %i' % nentries
        print 'ObsIDs = (%i,%i): %i' % (ObsIDextrema+(nObsID,))
        print 'Nr. Tests: %i' % ntests
        
        for it,test in enumerate(tests):
            print '%s\t %.2f min' % (test,durations[it])
        
        print 'Total Run Time = %.2f h' % totdurh
        
        


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
    
