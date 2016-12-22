#! /unsafe/raf/SOFTWARE/anaconda/envs/VISSIM/bin/python
# -*- coding: utf-8 -*-
"""

House-Keeping inspection tools.

:History:
Created on Thu Mar 10 12:11:58 2016

@author: Ruyman Azzollini
"""

# IMPORT STUFF
from pdb import set_trace as stop
import os
from glob import glob
import datetime

from matplotlib import pyplot as plt
import pylab
import matplotlib

from astropy.io import ascii
import string as st
import numpy as np
# END IMPORT


allHK_keys = {'5.7.02':['HK_OD_Top_CCD1','HK_OD_Bottom_CCD1','HK_OD_Top_CCD2',
'HK_OD_Bottom_CCD2','HK_OD_Top_CCD3','HK_OD_Bottom_CCD3','HK_RD_Top_CCD1',
'HK_RD_Bottom_CCD1','HK_RD_Top_CCD2','HK_RD_Bottom_CCD2','HK_RD_Top_CCD3',
'HK_RD_Bottom_CCD3','HK_temp_top_CCD1','HK_temp_bottom_CCD1','HK_temp_top_CCD2',
'HK_temp_bottom_CCD2','HK_temp_top_CCD3','HK_temp_bottom_CCD3','HK_IG1_top',
'HK_IG1_bot','HK_IG2_top','HK_IG2_bot','HK_IDH','HK_IDL','HK_DD_bias',
'HK_OG_bias','HK_1.5V_ROE','HK_VCCD_ROE','HK_5VA_pos_ROE','HK_5V_ref_ROE',
'HK_10VA_ROE','HK_5.2V_neg_ROE','HK_3V_neg_ROE','HK_VRclk_ROE',
'HK_VRClk_Lo_ROE','HK_3.3V_DIG_RPSU','HK_I3.3V_DIG_RPSU','HK_1.5V_DIG_RPSU',
'HK_I1.5V_DIG_RPSU','HK_28V_Pri_RPSU','HK_I28V_RPSU','HK_VAN_pos_RPSU',
'HK_I+VAN_RPSU','HK_VAN_neg_RPSU','HK_I-VAN_RPSU','HK_VCLK_RPSU',
'HK_IVCLK_RPSU','HK_VCCD_RPSU','HK_IVCCD_RPSU','HK_Temp1_RPSU',
'HK_3.3V_ROE','HK_Temp2_RPSU','HK_VClk_ROE'],
'5.7.04':['HK_OD_Top_CCD1','HK_OD_Bottom_CCD1','HK_OD_Top_CCD2',
'HK_OD_Bottom_CCD2','HK_OD_Top_CCD3','HK_OD_Bottom_CCD3','HK_RD_Top_CCD1',
'HK_RD_Bottom_CCD1','HK_RD_Top_CCD2','HK_RD_Bottom_CCD2','HK_RD_Top_CCD3',
'HK_RD_Bottom_CCD3','HK_temp_top_CCD1','HK_temp_bottom_CCD1','HK_temp_top_CCD2',
'HK_temp_bottom_CCD2','HK_temp_top_CCD3','HK_temp_bottom_CCD3','HK_IG1_top',
'HK_IG1_bot','HK_IG2_top','HK_IG2_bot','HK_IDH','HK_IDL','HK_DD_bias',
'HK_OG_bias','HK_1.5V_ROE','HK_VCCD_ROE','HK_5VA_pos_ROE','HK_5V_ref_ROE',
'HK_10VA_ROE','HK_5.2V_neg_ROE','HK_3V_neg_ROE','HK_VRclk_ROE',
'HK_VRClk_Lo_ROE','HK_3.3V_DIG_RPSU','HK_I3.3V_DIG_RPSU','HK_1.5V_DIG_RPSU',
'HK_I1.5V_DIG_RPSU','HK_28V_Pri_RPSU','HK_I28V_RPSU','HK_VAN_pos_RPSU',
'HK_I+VAN_RPSU','HK_VAN_neg_RPSU','HK_I-VAN_RPSU','HK_VCLK_RPSU',
'HK_IVCLK_RPSU','HK_VCCD_RPSU','HK_IVCCD_RPSU','HK_Temp1_RPSU',
'HK_3.3V_ROE','HK_Temp2_RPSU','HK_VClk_ROE'],
'5.7.07':['HK_OD_Top_CCD1','HK_OD_Bottom_CCD1','HK_OD_Top_CCD2',
'HK_OD_Bottom_CCD2','HK_OD_Top_CCD3','HK_OD_Bottom_CCD3','HK_RD_Top_CCD1',
'HK_RD_Bottom_CCD1','HK_RD_Top_CCD2','HK_RD_Bottom_CCD2','HK_RD_Top_CCD3',
'HK_RD_Bottom_CCD3','HK_temp_top_CCD1','HK_temp_bottom_CCD1','HK_temp_top_CCD2',
'HK_temp_bottom_CCD2','HK_temp_top_CCD3','HK_temp_bottom_CCD3','HK_IG1_top',
'HK_IG1_bot','HK_IG2_top','HK_IG2_bot','HK_IDH','HK_IDL','HK_DD_bias',
'HK_OG_bias','HK_1.5V_ROE','HK_VCCD_ROE','HK_5VA_pos_ROE','HK_5V_ref_ROE',
'HK_10VA_ROE','HK_5.2V_neg_ROE','HK_3V_neg_ROE','HK_VRclk_ROE','HK_VRClk_Lo_ROE',
'HK_3.3V_DIG_RPSU','HK_I3.3V_DIG_RPSU','HK_1.5V_DIG_RPSU','HK_I1.5V_DIG_RPSU',
'HK_28V_Pri_RPSU','HK_I28V_RPSU','HK_VAN_pos_RPSU','HK_I+VAN_RPSU',
'HK_VAN_neg_RPSU','HK_I-VAN_RPSU','HK_VCLK_RPSU','HK_IVCLK_RPSU','HK_VCCD_RPSU',
'HK_IVCCD_RPSU','HK_Temp1_RPSU','HK_3.3V_ROE','HK_Temp2_RPSU','HK_VClk_ROE',
'Video_TOP','Video_BOT']}

allHK_keys['5.7.06'] = allHK_keys['5.7.04']
allHK_keys['5.7.08'] = allHK_keys['5.7.07']
allHK_keys['5.7.09'] = allHK_keys['5.7.08']

allstats = ['mean','std','min','max']


def loadHK(filename,form='5.7.07'):
    """Loads a HK file
    
    It only assumes a structure given by a HK keyword followed by a number of 
    of blank-separated values (number not specified).
    Note that the length of the values arrays is variable (depends on length
    of exposure and HK sampling rate).
    
    :param filename: path to the file to be loaded, including the file itself
    :param form: format of HK file, given by version of "ELVIS"
    
    :return: dictionary with pairs parameter:[values]
    
    """
    f = open(filename)
    lines = f.readlines()
    f.close()
    
    data = {}
    
    for line in lines:
        items = st.split(line)
        key = items[0]
        values = np.array(items[1:]).astype('float32')
        
        data[key] = values
    
    
    return data

def filtervalues(values,key):
    """ """
    
    if ('HK_temp_top_CCD' in key) or ('HK_temp_bottom_CCD' in key):
        defvalue = -246.86
    elif ('HK_Temp1_RPSU' in key) or ('HK_Temp2_RPSU' in key):
        defvalue = -273.27
    elif (key == 'Video_TOP') or (key=='Video_BOT'):
        defvalue = -273.15
    else:
        defvalue = 0.
    
    try: return values[np.where(values != defvalue)]
    except: return np.zeros_like(values)+np.nan


def synthHK(HK):
    """ 
    Synthetizes the values for each parameter in a HK dictionary into
    [mean,std,min,max].
    
    
    :param HK: a dictionary as those output by loadHK.
    
    :return: dictionary with pairs parameter:[mean,std,min,max]
    
    """
    
    synthHKdict = {}
    
    for key in HK.keys():
        values = HK[key].copy()
        
        values = filtervalues(values,key)
        
        try: mean = values.mean()
        except: mean = np.nan
        try: std = values.std()
        except: std = np.nan 
        try: minv=  values.min()
        except: minv= np.nan
        try: maxv = values.max()
        except: maxv = np.nan
        
        synthHKdict[key] = np.array([mean,std,minv,maxv])

    return synthHKdict

def reportHK(HKs,key,reqstat='all'):
    """Returns (mean, std, min, max) for each keyword in a list of
    HK dictionaries (output from loadHK).
    
    :param HK: dictionary with HK data.
    :param key: HK key.
    :reqstat: what statistic to retrieve.
    
    """
    
    if reqstat != 'all': ixstat = allstats.index(reqstat)
    else: ixstat = np.where(allstats)
    
    vals = []
    
    for HK in HKs:
        vals.append(synthHK[HK][ixstat])
    
    vals = np.array(vals)
    
    return vals

def parseHKfname(HKfname):
    """Parses name of a HK file to retrieve OBSID, date and time, and ROE number.
    
    :parameter HKfname: name of HK file.
    :return: obsid,dtobj=datetime.datetime(yy,MM,dd,hh,mm,ss),ROE
    
    
    """
    
    stripHKfname = os.path.basename(HKfname)
    stripHKfname = os.path.splitext(stripHKfname)[0]
    
    items = st.split(stripHKfname,'_')
    
    assert items[0] == 'HK' # just a basic check
    
    obsid = int(items[1])
    DT = items[2]
    ROE = items[3]
    
    date = DT[0:DT.index('D')]
    dd,MM,yy = int(date[0:2]),int(date[2:4]),int(date[4:6])+2000

    time = DT[DT.index('D')+1:-1]
    hh,mm,ss = int(time[0:2]),int(time[2:4]),int(time[4:6])
    
    DTobject = datetime.datetime(yy,MM,dd,hh,mm,ss)
    
    return obsid,DTobject,ROE
    

def parseHKfiles(HKlist,form='5.7.07'):
    """ 
    
    :param HKlist: list of HK files (path+name).
    :return: [obsids],[dtobjs],[tdeltasec],[HK_keys], [data(nfiles,nstats,nHKparams)]
    
    """
    
    nfiles = len(HKlist)
    HK_keys = allHK_keys[form]
    nkeys = len(HK_keys)
    
    
    obsids = []
    dtobjs = []
    
    data = np.zeros((nfiles,4,nkeys),dtype='float32')
    
    for ix,HKfname in enumerate(HKlist):
        
        obsid,dtobj,ROE = parseHKfname(HKfname)        
        obsids.append(obsid)
        dtobjs.append(dtobj)
        
        HK = loadHK(HKfname)
        synthHKdict = synthHK(HK)
        
        for ik,key in enumerate(HK_keys):
            data[ix,:,ik] = synthHKdict[key]
        
    
    # reordering
    
    tsort = np.argsort(dtobjs)
    obsids = np.array(obsids)[tsort]
    dtobjs = np.array(dtobjs)[tsort]
    data = np.array(data)[tsort,:,:]
    
    tdeltasec = np.array([(item-dtobjs[0]).seconds for item in dtobjs])
    
   
    return obsids,dtobjs,tdeltasec,HK_keys, data


def format_date(x,pos=None):
    try: return pylab.num2date(x).strftime('%H:%M:%S')
    except: pass
    #except: stop()

def HKplot(allHKdata,keylist,key,dtobjs,filename='',stat='mean'):
    """Plots the values of a HK parameter as a function of time.
    
    :parameter allHKdata: HKdata = [(nfiles,nstats,nHKparams)]
    :parameter keylist: list with all HK keys.
    :parameter key: selected key.
    :parameter tdeltahour: time axis.
    
    :return: None!!
    
    """
    
    ixkey = keylist.index(key)
    
    ixstat = allstats.index(stat)
    
    HKvals = allHKdata[:,ixstat,ixkey]
    

    ylabel = 'V'
        
    curr_keys = ['hk_vrclk_lo_roe','hk_i3.3v_dig_rpsu','hk_i1.5v_dig_rpsu',
                 'hk_i28v_rpsu','hk_i+van_rpsu','hk_i-van_rpsu',
                 'hk_ivclk_rpsu','hk_ivccd_rpsu','hk_vclk_roe']
   
    if key.lower() in curr_keys:
        ylabel = 'I[A]'
        
    temp_keys = ['hk_temp_top_ccd1','hk_temp_bottom_ccd1',
                 'hk_temp_top_ccd2','hk_temp_bottom_ccd2',
                 'hk_temp_top_ccd3','hk_temp_bottom_ccd3',
                 'hk_temp1_rpsu','hk_temp2_rpsu',
                 'video_top','video_bot']
    
    if key.lower() in temp_keys:
        ylabel = 'deg [C]'    

    fig = plt.figure(figsize=(7,6))
    ax1 = fig.add_subplot(111)
    if not np.all(np.isnan(HKvals)):
        ax1.plot(dtobjs,HKvals,'b-')
    #ax1.plot(HKvals,'b')
    #ax1.set_xlabel(r'$Time$')
    ax1.set_ylabel(r'$%s$' % ylabel)
    ax1.set_title(key)
    
    for item in [ax1.title, ax1.xaxis.label, ax1.yaxis.label]:
        item.set_fontsize(18)
    
    ylim = ax1.get_ylim()
    dy = ylim[1]-ylim[0]
    if dy < 0.1:
        ymid = (ylim[1]+ylim[0])/2.
        ylim = (ymid-0.05,ymid+0.05)
    
    ax1.set_ylim(ylim)
    
    #if key == 'HK_temp_top_CCD1': stop()
        
    #if dtobjs[0] != 0:
    plt.gca().xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(format_date)) 
        
    fig.autofmt_xdate()
    
    
    plt.tight_layout(rect=[0, 0, 1, 1])
    
    if filename== '':
        plt.show()
    else:
        plt.savefig(filename)
   
        
    plt.close()

def test():
    
    datapath = 'data/08_Mar/'
    HKfile = os.path.join(datapath,'HK_281_080316D174300T_ROE1.txt')
    
    #HK = loadHK(HKfile)
    #HKbrief = synthHK(HK)
    
    HKlist = glob('%s/HK_*ROE1.txt' % datapath)
    
    obsids, dtobjs, tdeltasec, HK_keys, allHKdata = parseHKfiles(HKlist)
    # [allHKdata] = [n_obsids, n_stats, n_hk]
    tdeltahour = tdeltasec / 3600.
    HK_temp_top_CCD2 = allHKdata[:,0,HK_keys.index('HK_temp_top_CCD2')]
    HK_temp_bot_CCD2 = allHKdata[:,0,HK_keys.index('HK_temp_bottom_CCD2')]
    
    
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax1.plot(tdeltahour,HK_temp_bot_CCD2,'b-',label='bottom')
    ax1.plot(tdeltahour,HK_temp_top_CCD2,'r-',label='top')
    ax1.set_xlabel(r'$\Delta T[hour]$')
    ax1.set_ylabel('Temps CCD2 [C]')
    ax1.set_ylim([-120.,-40.])
    
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles,labels)
    
    ax2 = fig.add_subplot(122)
    ax2.plot(tdeltahour,HK_temp_bot_CCD2,'b-')
    ax2.plot(tdeltahour,HK_temp_top_CCD2,'r-')
    ax2.set_xlabel(r'$\Delta T[hour]$')
    ax2.set_ylabel('Temps CCD2 [C]')
    #ax2.set_ylim([-105,-95])
    
    fig.suptitle('08.March.2016 - Functionals')
    plt.show()
    
    
    
    for key in HK_keys:
        HKplot(allHKdata,HK_keys,key,tdeltahour)
    
if __name__ == '__main__': 
    
    test()
