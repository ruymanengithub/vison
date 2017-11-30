# -*- coding: utf-8 -*-
"""
Script to analyze Test PSF0X

PSF vs. Fluence, and Wavelength
   PSF01 - nominal temperature
   PSF02 - alternative temperatures

Tasks:

    - select exposures, get file names, extract data.
    - subtract offset level.
    - divide by Flat-field.
    - crop stamps of the sources on each CCD/Quadrant.
       - save snapshot figures of sources.
    - for each source:
       - measure shape using weighted moments
       - measure shape using Gaussian Fit
       - Forward Model the optomechanic+detector PSF

Created on Fri Nov 25 19:14:12 2016

@author: raf
"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
from collections import OrderedDict

from vison.plot import classes as plclasses

# END IMPORT

class genplot_stats_beam:
    """ """



PSF0Xfigs = dict()
PSF0Xfigs['P0Xchecks_offsets'] = genplot_stats_beam
PSF0Xfigs['P0Xchecks_stds'] = genplot_stats_beam
PSF0Xfigs['BlueScreen'] = plclasses.BlueScreen



#def subbgd(stamp):
#    """ """
#    
#    if isinstance(stamp,np.ma.masked_array):
#        stamp = stamp.data
#    
#    bgd = np.nanmedian(stats.sigma_clip(stamp,sigma=6))
#    stamp_mbgd = stamp - bgd
#    
#    return stamp_mbgd
#
#
#def showstamp(stamp,title='',outfile=None):
#    """ """
#    
#    fig = plt.figure()
#    ax = fig.add_subplot(111)
#    vmin = np.nanmedian(stamp)
#    vmax = np.max(stamp)
#    cax = ax.imshow(stamp.transpose(),origin='lower left',
#                    vmin=vmin,vmax=vmax)
#    cax.set_cmap('spectral')
#    ax.set_title(title)
#    cbar = fig.colorbar(cax,orientation='vertical')
#    if outfile is not None:
#        plt.savefig(outfile)
#    else:
#        plt.show()
#    plt.close()
#    
#def process_exposures_PSF0X(DataDict,box,Quad,mask_fits,flatfields_dict,VSCAN=[500-1,800-1],doshow=False):
#    """ 
#    Flow:
#    
#
#     - load (wavelength-matching) FF
#     
#     - for each exposure:
#          subtract offset level (from overscan)
#          divide by FF
#          cut stamp
#          save stamp with its metadata
#          
#     - return updated DataDict
#    
#    """
#    
#    wavelength = DataDict['wavelength']
#    wavekey = 'nm%i' % wavelength
#    
#    # Windowing
#    
#    dumbccd = ccdmod.CCD()
#    VscanMask = dumbccd.do_Vscan_Mask(VSCAN[0],VSCAN[1]).copy()
#    
#    
#    # LOAD MASK
#    
#    mask = fts.getdata(mask_fits).transpose().astype('bool')
#    
#    mask = mask | VscanMask    
#
#    
#    # LOAD FF
#    
#
#    FF_fits = flatfields_dict[wavekey]
#    FF = fts.getdata(FF_fits).transpose() # [:,:,0]
#
#    # Exposure by Exposure Processing
#        
#    uexptime = DataDict['uexptime']
#    
#    expkeys = ['ms%i' % item for item in uexptime]
#    
#    for iexp,expkey in enumerate(expkeys):
#        
#        print '%s-%s... %i/%i exposure times' % (wavekey,expkey,iexp+1,len(expkeys))
#        
#        Files = DataDict[expkey]['Files']
#        OBSIDS = DataDict[expkey]['OBSIDS']
#        Exptime = DataDict[expkey]['Exptime']
#        
#        DataDict[expkey]['SHUTTER'] = np.zeros_like(Files,dtype='bool')
#        DataDict[expkey]['PICKLES'] = np.zeros_like(Files,dtype='40str')
#        
#        
#        for iobs,iinfits in enumerate(Files):
#            
#            print 'Exposure %i / %i' % (iobs+1,len(Files))
#            
#            #if OBSIDS[iobs] != 2591: continue
#            
#            OBSID = OBSIDS[iobs]  
#            exptime = Exptime[iobs]
#            iinfits = os.path.join(datapath,iinfits)
#    
#            imgobj = ccdmod.CCD(iinfits)
#            
#            # Masking
#            imgobj.get_mask(mask)
#            
#            
#            # Offset subtraction
#            
#            for Q in Quads:
#                imgobj.sub_offset(Q,method='median',scan='ove',trimscan=[4,3])
#            
#            
#            # FF Division
#
#            imgobj.data /= FF
#            
#    
#            imgdata = imgobj.get_quad(Quad,canonical=False)
#            stamp = imgdata[box[0]:box[1],box[2]:box[3]]
#            if doshow: showstamp(stamp,title='%i:%s' % (OBSID,Quad))
#            
#            # BACKGROUND SUBTRACTION
#    
#            stamp_nobgd = subbgd(stamp) 
#            
#            # SAVING
#    
#            outfile = 'PSF0X_%i_Q%s_CCD%i_ofsflfbgd.cpickle' % (OBSID,Quad,CCD)
#            outfile = os.path.join(resultspath,outfile)
#            
#            FAILED = stamp_nobgd.max() < 1.E2 # failed exposure, no source
#
#            results = dict(img=stamp_nobgd,box=box,\
#              FF=os.path.split(flatfields_dict[wavekey])[-1],
#              MASK=os.path.split(mask_fits)[-1],
#              RON=FMlib.RON,gain=FMlib.gain,exptime=exptime,FAILED=FAILED) 
#    
#            cPickleDumpDictionary(results,outfile)
#            
#            DataDict[expkey]['SHUTTER'][iobs] = ~FAILED
#            DataDict[expkey]['PICKLES'][iobs] = os.path.split(outfile)[-1]
#    
#    return DataDict

#def filter_exposuresPSF01(explog,OBSID_lims,inwavelength,inCCD):
#    """ """
#    
#    OBSIDs = np.array(explog['ObsID'].copy())
#    #Wavelength = np.array(explog['Wavelength'].copy()) # EXP_LOG corrupted, wavelenth selection not possible
#    Exptime = np.array(explog['Exptime'].copy())
#    rootFile_name = np.array(explog['File_name'].copy())
#    TEST = np.array(explog['TEST'].copy())
#    CCD = np.array(explog['CCD'].copy())
#    
#    DataDict = {}
#        
#    DataDict = dict(wavelength=inwavelength)
#    #selbool = (TEST == 'PSF02') & (Wavelength == inwavelength) & \ # IF THE LOG WAS RIGHT...
#    #    (CCD == 'CCD%i' % inCCD)
#
#    selbool = (TEST == 'PSF01') & (OBSIDs >= OBSID_lims[0]) & \
#        (OBSIDs <= OBSID_lims[1]) & (CCD == 'CCD%i' % inCCD) & (Exptime > 0.)
#    
#    
#    ixsel = np.where(selbool)
#    uexptime = np.unique(Exptime[ixsel])
#    
#    nlevels = len(uexptime)
#    DataDict['uexptime'] = uexptime
#    
#    
#    for iexp in range(nlevels):
#        
#        expkey = 'ms%i' % uexptime[iexp]
#        DataDict[expkey] = {}
#        
#        ixsel = np.where(selbool & (Exptime == uexptime[iexp]))
#        
#        DataDict[expkey]['OBSIDS'] = OBSIDs[ixsel]
#        
#        ixrootFile_name = rootFile_name[ixsel]
#        
#        ixFile_name = ['%s.fits' % item for item in ixrootFile_name]
#        
#        DataDict[expkey]['Files'] = ixFile_name
#        DataDict[expkey]['Exptime'] = Exptime[ixsel]
#        
#    
#    return DataDict



#def run(inputs,log=None):
#    
#    
#    if log: log.info('PSF01 Running...')
#    
#    
#    sys.exit('PSF01 not ready')
#    
#    FilterID = inputs['FilterID']
#    VSTART = inputs['VSTART']
#    VEND = inputs['VEND']
#    doshow =inputs['doshow']
#    box = inputs['box']
#    masks_dict = inputs['masks_dict']
#    flatfields_dict = inputs['flatfields_dict']
#    explogf = inputs['explogf']
#    resultsroot = inputs['resultsroot']
#    datapath = inputs['datapath']
#    
#    resultspath = os.path.join(resultsroot,'PSF01')
#    
#    global resultspath
#    global datapath
#    
#    
#    wavelength = FMlib.FW[FilterID]
#    wavekey = '%inm' % wavelength
#    
#    outdict = os.path.join(resultspath,'DataDict_PSF01_%s.cpickle' % wavekey)
#    
#    explog = ELtools.loadExpLog(explogf)
#    
#    OBSID_lims = [OBSID_start,OBSID_end]
#    
#    DataDict = filter_exposuresPSF01(explog,OBSID_lims,wavelength,CCD)
#    
#    DataDict = process_exposures_PSF01(DataDict,box,Q,mask_fits,\
#       flatfields_dict=flatfields_dict,VSCAN=[VSTART,VEND],doshow=doshow)
#    
#    #wavekeys = DataDict.keys()
#    #for wavekey in wavekeys:
#        
#    uexptime = DataDict['uexptime']
#        
#    expkeys = ['ms%i' % item for item in uexptime]
#        
#    for expkey in expkeys:
#            
#        SHUTTER = DataDict[expkey]['SHUTTER'].copy()
#        Nfailed = len(np.where(SHUTTER == False)[0])
#            
#        print '%s-%s: %i exposures FAILED out of %i' % \
#             (wavekey,expkey,Nfailed,len(SHUTTER))
#    
#    cPickleDumpDictionary(DataDict,outdict)
#
#    return None
    
    
#if __name__ == '__main__':
#    
#    doRun = False
#    wavelengths = [570,700,800,890]
#    #flevels = [1,2,3,4,5]
#    VSTART = 500-1
#    VEND = 800-1
#    doshowstamps = False
#    box = [1103,1201,600,686]
#    Q= 'H'
#    CCD = 2
#    
#    OBSIDs = dict(nm570=[2099,2197],nm700=[1985,2083],nm800=[1864,1962],nm890=[1740,1838])
#        
#    
#    mask_fits = os.path.join('COSMETICS','results','Defects_CharEM1A_sn14173-11-02.fits')
#    #superbias_fits = 'SuperBias_CCD273.fits'
#    flatfields_dict = dict(nm570='EUC_FFw545nm_IL2.fits',\
#                           nm700='EUC_FFw750nm_IL2.fits',nm800='EUC_FFw800nm_IL2.fits',\
#                           nm890='EUC_FFw850nm_IL2.fits')
#    for fkey in flatfields_dict.keys():
#        flatfields_dict[fkey] = os.path.join('BEST_FLATFIELDS','results_spline',flatfields_dict[fkey])
#    
#    
#    if doRun:
#    
#        inputs = {}
#        inputs['VSTART'] = VSTART
#        inputs['VEND'] = VEND
#        inputs['doshow'] = doshowstamps
#        inputs['box'] = box
#        inputs['Q'] = Q
#        inputs['CCD'] = CCD
#        inputs['mask_fits'] = mask_fits
#        inputs['flatfields_dict'] = flatfields_dict
#        
#        for wavelength in wavelengths:
#            
#            wavekey = 'nm%i' % wavelength
#            
#            if wavelength != 570:
#                datapath = '21_Oct_16'
#                datakey = '211016'
#            else:
#                datapath = '23_Oct_16'
#                datakey = '231016'
#                    
#            resultspath = 'results_nm%i' % wavelength
#            if not os.path.exists(resultspath):
#                os.system('mkdir %s' % resultspath)
#            
#        
#            explogf = os.path.join(datapath,'EXP_LOG_%s.txt' % datakey)
#            
#            inputs['explogf'] = explogf
#            inputs['OBSID_start'] = OBSIDs[wavekey][0]
#            inputs['OBSID_end'] = OBSIDs[wavekey][1]
#            
#            print '\nProcessing Wavelength = %i nm\n' % wavelength
#            
#            run(wavelength,inputs)
#            
#    else:
#        
#        for wavelength in wavelengths:
#            
#            wavekey = 'nm%i' % wavelength
#            
#            if wavelength != 570:
#                datapath = '21_Oct_16'
#                datakey = '211016'
#            else:
#                datapath = '23_Oct_16'
#                datakey = '231016'
#                    
#            resultspath = 'results_nm%i' % wavelength
#            if not os.path.exists(resultspath):
#                os.system('mkdir %s' % resultspath) 
#        
#            wavekey = '%inm' % wavelength
#    
#            outdict = os.path.join(resultspath,'DataDict_PSF01_%s.cpickle' % wavekey)
#            
#            DataDict = cPickleRead(outdict)
#            uexptime = DataDict['uexptime']
#            
#            for utime in uexptime:
#                utimekey = 'ms%i' % utime
#                OBSIDS = DataDict[utimekey]['OBSIDS']
#                PICKLES = DataDict[utimekey]['PICKLES']
#                
#                nOBS = len (DataDict[utimekey]['OBSIDS'])
#                
#                for ix in range(nOBS):
#                    
#                    picklef = os.path.join(resultspath,PICKLES[ix])
#                    
#                    pickle = cPickleRead(picklef)
#                    
#                    stamp = pickle['img']
#                    
#                    stampf = os.path.join(resultspath,'Stamp_%i_ms%i_nm%i.png' % \
#                        (OBSIDS[ix],utime,wavelength))
#                    
#                    title = 'OBSID-%i, ms%i, nm%i' % (OBSIDS[ix],utime,wavelength)                    
#                    
#                    showstamp(stamp,title,stampf)
#            
            