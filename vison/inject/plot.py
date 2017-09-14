#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Charge Injection Plotting Tools.

Created on Thu Sep 14 15:39:34 2017

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk


"""

# IMPORT STUFF

from matplotlib import pyplot as plt

import numpy as np
from pdb import set_trace as stop

import matplotlib
matplotlib.rcParams['font.size'] = 17
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('axes', linewidth=1.1)
matplotlib.rcParams['legend.fontsize'] = 12
matplotlib.rcParams['legend.handlelength'] = 3
matplotlib.rcParams['xtick.major.size'] = 5
matplotlib.rcParams['ytick.major.size'] = 5
matplotlib.rcParams['image.interpolation'] = 'none'
import matplotlib.cm as cm

# END IMPORT


def get_phaseV(phaseID,delay,TOI=500):
    """A1=Iphi2
    A = [G,H]
    TOI=500 us
    """
    
    waveform = [0,-3,8,8]
    
    steps0 = dict(A1=3,A2=2,A3=1,A4=0,
                D1=0,D2=3,D3=2,D4=1)
    
    
    step0 = steps0[phaseID]
    
    step = ((delay / TOI) + step0) % 4
    
    V = waveform[step]
    return V
    
    

def draw_gate_voltages(IG1,IG2,IDL,IDH,delay,figname='',chpar=10,rchpar=7.5,TOI=500):
    """ 
    IG1: float, volts
    IG2: float, volts
    IDL: float, volts
    IDH: float, volts
    delay : float, microseconds
    
    """
    

    upperphases = ['A1','A2','A3','A4'] # E&F
    lowerphases = ['D1','D2','D3','D4'] # G&H

    unity = np.ones(2)
    stepx = np.array([-0.5,+0.5])
    halfstepx = np.array([0.,+0.5])
    
    
    gatenames = ['IG1','IG2']
    
    phaselabels = dict(D1=r'$D1=I\phi1$',D2=r'$D2=I\phi2$',D3=r'$D3=I\phi3$',D4=r'$D4=I\phi4$',
                       A1=r'$A1=I\phi2$',A2=r'$A2=I\phi3$',A3=r'$A3=I\phi4$',A4=r'$A4=I\phi1$')
    
    fig = plt.figure()
    
    yrange = [5,23]
    
    static_voltages = [IG1,IG2]
    
    ax1 = fig.add_subplot(211)
       
    ax1.plot(0+stepx,IDH*unity,'-',c='b',lw=2)
    ax1.text(0,IDH+0.5,'ID',fontsize=12,verticalalignment='top')
    
    #IDLcolors = cm.rainbow( np.linspace( 0,1,len(IDL)))
    
    #for ix,iIDL in enumerate(IDL):
    #    ax1.plot(0+stepx,iIDL*unity,'--',c=IDLcolors[ix],lw=2,label='IDL=%.1fV' % iIDL)
    
    # IDL & IDH
    
    ax1.plot(0+stepx,IDL*unity,'--',c='r',lw=2)
    ax1.plot(0+stepx,IDH*unity,'-',c='b',lw=2)
    
    for iVgate,Vgate in enumerate(static_voltages):
        VV = rchpar+static_voltages[iVgate]
        
        if iVgate == 1:
            VVnotch = VV + chpar-rchpar
            ax1.plot(1+iVgate+halfstepx,VVnotch*unity,'-',c='g',lw=2)
            ax1.plot(1+iVgate-halfstepx,VV*unity,'-',c='g',lw=2)
            
            ax1.text(1+iVgate,VVnotch+0.5,gatenames[iVgate],fontsize=12,\
            verticalalignment='top')
        else:
            ax1.plot(1+iVgate+stepx,VV*unity,'-',c='g',lw=2)
            ax1.text(1+iVgate,VV+0.5,gatenames[iVgate],fontsize=12,\
            verticalalignment='top')
    
    for iVphase,Vphase in enumerate(upperphases):
        VV=chpar+get_phaseV(Vphase,delay,TOI=TOI)
        ax1.plot(3+iVphase+stepx,VV*unity,'-',c='r',lw=2)
        
        ax1.text(3+iVphase-0.5,VV+0.5,phaselabels[Vphase],fontsize=12,\
        verticalalignment='top',horizontalalignment='left')
        
    
    for ix in np.arange(-0.5,7.5,1.):
        ax1.axvline(x=ix,color='k',ls='--')
    ax1.axvline(x=2,color='k',ls=':')
    
    ax1.axes.get_xaxis().set_visible(False)
    ax1.set_title('E & F')
    
    ax1.set_ylim(yrange)
    ax1.set_ylim(ax1.get_ylim()[::-1])
    ax1.set_ylabel('Eff. Voltage')
    
    
    ax2 = fig.add_subplot(212)
       
    ax2.plot(0+stepx,IDH*unity,'-',c='b',lw=2)
    ax2.text(0,IDH+0.5,'ID',fontsize=12,verticalalignment='top')
    
    # IDL & IDH
    
    ax2.plot(0+stepx,IDL*unity,'--',c='r',lw=2)
    ax2.plot(0+stepx,IDH*unity,'-',c='b',lw=2)
    
    #for ix,iIDL in enumerate(IDL):
    #    ax2.plot(0+stepx,iIDL*unity,'--',c=IDLcolors[ix],lw=2)   
    #ax2.plot(0+stepx,IDH*unity,'-',c='b',lw=2)
    
    for iVgate,Vgate in enumerate(static_voltages):
        VV = rchpar+static_voltages[iVgate]
        
        if iVgate == 1:
            VVnotch = VV + chpar-rchpar
            ax2.plot(1+iVgate+halfstepx,VVnotch*unity,'-',c='g',lw=2)
            ax2.plot(1+iVgate-halfstepx,VV*unity,'-',c='g',lw=2)
            ax2.text(1+iVgate,VVnotch+0.5,gatenames[iVgate],fontsize=12,\
            verticalalignment='top')
        else:
            ax2.plot(1+iVgate+stepx,VV*unity,'-',c='g',lw=2)
            ax2.text(1+iVgate,VV+0.5,gatenames[iVgate],fontsize=12,\
            verticalalignment='top')
    
    
    for iVphase,Vphase in enumerate(lowerphases):
        VV=chpar+get_phaseV(Vphase,delay)
        
        ax2.plot(3+iVphase+stepx,VV*unity,'-',c='r',lw=2)
        
        ax2.text(3+iVphase-0.5,VV+0.5,phaselabels[Vphase],fontsize=12,\
        verticalalignment='top',horizontalalignment='left')
    
    for ix in np.arange(-0.5,7.5,1.):
        ax2.axvline(x=ix,color='k',ls='--')
    ax2.axvline(x=2,color='k',ls=':')
    
    
    ax2.axes.get_xaxis().set_visible(False)
    ax2.set_title('G & H')
    
    ax2.set_ylim(yrange)
    ax2.set_ylim(ax2.get_ylim()[::-1])
    ax2.set_ylabel('Eff. Voltage')
    
    title = 'IDL=%.1fV IDH=%.1fV IG1=%.1fV IG2=%.1fV Del=%i us' % \
     (IDL,IDH,IG1,IG2,delay)
    
    plt.suptitle(title)
    
    #lines,handles = ax1.get_legend_handles_labels()
    #fig.legend(lines,handles,'right')
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.85)
        
    
    if figname != '':
        plt.savefig(figname)
    else:
        plt.show()
    plt.close()