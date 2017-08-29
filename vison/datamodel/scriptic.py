#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Classes and functions to generate ELVIS commanding scripts automatically.

Created on Wed May 24 15:31:54 2017

:author: Ruyman Azzollini
:contact: r.azzollini__at__ucl.ac.uk

"""

# IMPORT STUFF
from pdb import set_trace as stop
import pandas as pd
import string
import numpy as np
# END IMPORT


script_dictionary = {'6.0.0':{'keys':['frames',
'program','test','IDL','IDH','IG1','IG2','OD_1',
'RD_1','OD_2','RD_2','OD_3','RD_3',
'iphi1','iphi2','iphi3','iphi4',
'readmode_1','readmode_2','vertical_clk','serial_clk',
'flushes','exptime','shutter','electroshutter','vstart','vend',
'sinvflush',
'chinj','chinj_rows_on','chinj_rows_off','chinj_repeat',
'id_width','id_delay',
'tpump','ser_shuffles','ver_shuffles',
'dwell_v','dwell_h',
'motor','matrix_size','step_size',
'add_h_overscan','add_v_overscan',
'toi_flush','toi_tpump','toi_rdout','toi_chinj',
'wavelength','pos_cal_mirror',
'operator',
'sn_ccd1','sn_ccd2','sn_ccd3','sn_roe','sn_rpsu',
'comments'],
'alias':dict(frames='Frames',
program='Program',test='Test',IDL='IDL(mV)',IDH='IDH(mV)',
 IG1='IG1(mV)',IG2='IG2(mV)',
OD_1='OD CCD1(mV)',RD_1='RD CCD1(mV)',OD_2='OD CCD2(mV)',
             RD_2='RD CCD2(mV)',OD_3='OD CCD3(mV)',RD_3='RD CCD3(mV)',
iphi1='Iphi1',iphi2='Iphi2',iphi3='Iphi3',iphi4='Iphi4',
readmode_1='Readout mode R01',readmode_2='Readout mode R02',
vertical_clk='Vertical clk',serial_clk='Serial clk',
flushes='Flushes',exptime='Exposure time(ms)',
shutter='Shutter',electroshutter='Electronic shutter',
vstart='Vstart',vend='Vend',
sinvflush='S. Inv.Flush',
chinj='charge injection',chinj_rows_on='chrg_inj rows on',
chinj_rows_off='chrg_inj rows off',chinj_repeat='chrg_inj repeat',
id_width='id pulse length(us)',id_delay='id pulse delay(us)',
tpump='Trap pumping',ser_shuffles='Serial shuffles',ver_shuffles='Vertical shuffles',
dwell_v='Dwell_V',dwell_h='Dwell_H',
motor='Move motor',matrix_size='Matrix size',step_size='Step size',
add_h_overscan='Additional H.Overscan',add_v_overscan='Additional V.Overscan',
toi_flush='TOI Flush(us)',toi_tpump='TOI T.Pump(us)',
                    toi_rdout='TOI Rdout(us)',toi_chinj='TOI Chinj(us)',
wavelength='Wavelength',pos_cal_mirror='Pos cal mirro(mm)',
operator='Operator name',
sn_ccd1='CCD1 type/sn',sn_ccd2='CCD2 type/sn',sn_ccd3='CCD3 type/sn',
sn_roe='ROE type/sn',sn_rpsu='RPSU type/sn',
comments='Comments')
    }}


class Script(object):
    """Core Class that provides automatic test script generation and validation."""
    
    def __init__(self,defaults={},structure={},elvis='6.0.0'):
        """Initialization."""
        
        self.scriptname = ''
        self.defaults = defaults
        self.structure = structure
        self.elvis = elvis
        self.cargo = []
    
    
    def build_cargo(self):
        """Updates 'cargo' attribute.
        'cargo': list of lists, each corresponding to a column in the script.
                 Each element in the inner lists is a register value.
                 The first column corresponds to the aliases column.
            
        """
        dictio = script_dictionary[self.elvis]
        keys = dictio['keys']
        aliases = dictio['alias']
        
        stru = self.structure
        defaults = self.defaults
        
        self.cargo = []
        
        self.cargo.append([aliases[key] for key in keys]) # first column, aliases
        
        cols = ['col%i' % (i+1,) for i in range(stru['Ncols'])]

        for col in cols:
            
            vallist = []
            
            for key in keys:
                if key in stru[col]:
                    val = stru[col][key]
                else:
                    val = defaults[key]
                vallist.append(val)
            
            self.cargo.append(vallist)
        
        if self.elvis > '6.0.0':
            self.cargo.append(['End'])
        
    
    def write(self,scriptname):
        """Writes self.cargo (script) to an excel file."""
        
        alphabet = list(string.ascii_lowercase)
        
        ncols = len(self.cargo)
        
        nalpha = ncols/len(alphabet)
        
        colnames = []
    
        colnameroot =''
        for ia in range(nalpha+1):
        
            if ia>0: colnameroot += alphabet[(ia-1)%len(alphabet)]
        
            for iia in range(len(alphabet)):
                colnames.append(colnameroot+alphabet[iia])
                if len(colnames) == ncols:
                    break
        
        datadict = {} 
        
        for ixc in range(ncols):
            datadict[colnames[ixc]] = self.cargo[ixc]
        
        df = pd.DataFrame(datadict)
 
        writer = pd.ExcelWriter(scriptname)
        df.to_excel(writer,'Sheet1',index=False)
        #writer.save()
        df.save()
        
        self.scriptname = scriptname

        return None
    
    def load(self,scriptname,elvis='6.0.0'):
        """Loads an script from an excel file."""
        self.elvis = elvis
        
        converters = {'Frames':str}
        df= pd.read_excel(scriptname,header=0,converters=converters)
        
        cargo = []
        cargo.append(['Frames']+[string.strip(item) for item in df['Frames'].tolist()])
        
        sdict = script_dictionary[self.elvis]
        keys = sdict['keys']
        alias = sdict['alias']
        aliaslist = [alias[key] for key in keys]
        

        assert np.all(cargo[0] == aliaslist)
        
        for ixc in df.columns:
            if ixc == 'Frames':
                continue
            elif ixc == 'End':
                break
            elif isinstance(ixc,int):
                rawvals = [ixc] + df[ixc].tolist()
                vals=[]
                for item in rawvals:
                    if isinstance(item,unicode):
                        vals.append(string.strip(str(item)))
                    else:
                        vals.append(item)
                
                cargo.append(vals)
    
        self.cargo = cargo
        
        self.scriptname = scriptname
        
        
    
    
    def validate(self,defaults,structure,elvis='6.0.0'):
        """Not sure 'validation' will work like as implemented...
        TODO: validate self.validate"""
        
        assert self.elvis == elvis
        
        access = Script(defaults,structure,elvis=elvis)
        access.build_cargo()
        
        assert len(access.cargo) == len(self.cargo)
        
        for i in range(len(self.cargo)):
            assert np.all(self.cargo[i] == access.cargo[i])
        
        
    
    

def test0():
    """ """
    
    scriptname = '../../../../CALIBRATION/CAMPAIGN/ELVIS/ELVIS_6.0.0_27Jan17/EUCLID_QM_ELVIS_script_v6.0.0.xlsx'
    
    script = Script()
    script.load(scriptname)
    stop()

if __name__ == '__main__':
    
    
    test0()