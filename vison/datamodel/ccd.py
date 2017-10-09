# -*- coding: utf-8 -*-
"""

Data model for Euclid-VIS CCDs (ground testing at MSSL)

Created on Fri Nov 13 17:42:36 2015

:Author: Ruyman Azzollini
"""

# IMPORT STUFF
from astropy.io import fits as fts
import numpy as np
import os
from pdb import set_trace as stop
import sys
# END IMPORT

isthere = os.path.exists

NAXIS1 = 4238
NrowsCCD = 2066
NAXIS2 = (NrowsCCD+20)*2 # 4132
prescan = 51
overscan = 20
#imgarea = [2048,2066]
RON = 1.4
gain = 3.1 # e/adu

QuadBound = dict(E=[0,NAXIS1/2,NAXIS2/2,NAXIS2],\
               F=[NAXIS1/2,NAXIS1,NAXIS2/2,NAXIS2],\
               G=[NAXIS1/2,NAXIS1,0,NAXIS2/2],\
               H=[0,NAXIS1/2,0,NAXIS2/2])

Quads = ['E','F','G','H']


HeadKeys = {'6.0.0':['EXTNAME','BUNIT','PROGRAM','OBJECT','OBSID','OPERATOR',
'FULLPATH','LAB_VER','CON_FILE','DATE','EXPTIME','FL_RDOUT','CI_RDOUT',
'N_P_HIGH','CHRG_INJ','ON_CYCLE''OFF_CYCL''RPEAT_CY''PLS_WDTH','PLS_DEL',
'SERRDDEL','TRAPPUMP','TP_SER_S','TP_VER_S','TP_DW_V','TP_DW_H','TOI_FLSH',
'TOI_PUMP','TOI_READ','TOI_CINJ','INVFLSHP','INVFLUSH','FLUSHES','VSTART',
'VEND','OVRSCN_H','CLK_ROE','CNVSTART','SUMWELL','INISWEEP','SPW_CLK','FPGA_VER',
'EGSE_VER','M_STEPS','M_ST_SZE','WAVELENG','MIRR_POS','CHMB_PRE','CCD1_SN',
'CCD2_SN','CCD3_SN','ROE_SN','CALSCRPT','COMMENTS','TMPCCD1T','TMPCCD1B',
'TMPCCD2T','TMPCCD2B','TMPCCD3T','TMPCCD3B','IDL_V','IDH_V','IG1_T1_V','IG1_T2_V',
'IG1_T3_V','IG1_B1_V','IG1_B2_V','IG1_B3_V','IG2_T_V','IG2_B_V','OD_T1_V',
'OD_T2_V','OD_T3_V','OD_B1_V','OD_B2_V','OD_B3_V','RD_T_V','RD_B_V']}
HeadKeys['6.3.0'] = [
'XTENSION','BITPIX',
'NAXIS','NAXIS1','NAXIS2',
'PCOUNT','GCOUNT',
'BZERO','BSCALE',
'EXTNAME','BUNIT',
'PROGRAM','OBJECT','OBSID','OPERATOR','FULLPATH','LAB_VER','CON_FILE',
'DATE','EXPTIME','FL_RDOUT','CI_RDOUT',
'IPHI','CHINJ','CHINJ_ON','CHINJ_OF','ID_WID','ID_DLY','CHIN_DLY',
'V_TPUMP','S_TPUMP','V_TP_MOD','S_TP_MOD','V_TP_CNT','S_TP_CNT',
'DWELL_V','DWELL_S',
'TOI_FL','TOI_TP','TOI_RO','TOI_CH',
'SIFLSH','SIFLSH_P','FLUSHES',
'VSTART','VEND','RDMODE','SWELLW','SWELLDLY','INISWEEP',
'SPW_CLK','FPGA_VER','EGSE_VER','MOTR_ON','MOTR_CNT','MOTR_SIZ',
'SOURCE','LAMBDA','MIRR_ON','MIRR_POS',
'SN_CCD1','SN_CCD2','SN_CCD3','SN_ROE','SN_RPSU',
'CALSCRPT',
'COMMENTS','R1C1_TT','R1C1_TB','R1C2_TT','R1C2_TB','R1C3_TT','R1C3_TT',
'IDL','IDH',
'IG1_1_T','IG1_2_T','IG1_3_T','IG1_1_B','IG1_2_B','IG1_3_B',
'IG2_T','IG2_B',
'OD_1_T','OD_2_T','OD_3_T','OD_1_B','OD_2_B','OD_3_B','RD_T','RD_B']

#Head_2_Explog_Dict = {'6.0.0':{'EXTNAME':None,'BUNIT':None,'PROGRAM':'PROGRAM',
#'OBJECT':'TEST','OBSID':None,'OPERATOR':'Operator',
#'FULLPATH':None,'LAB_VER':'Lab_ver','CON_FILE':'Con_file',
#'DATE':'DATE','EXPTIME':'Exptime','FL_RDOUT':'Flsh-Rdout_e_time',
#'CI_RDOUT':'C.Inj-Rdout_e_time','N_P_HIGH':'N_P_high',
#'CHRG_INJ':'Chrg_inj','ON_CYCLE':'On_cycle','OFF_CYCL':'Off_cycle',
#'RPEAT_CY':'Rpeat_cy','PLS_WDTH':'pls_len','PLS_DEL':'pls_del',
#'SERRDDEL':'SerRDel','TRAPPUMP':'Trappump',
#'TP_SER_S':'TP_Ser_S','TP_VER_S':'TP_Ver_S','TP_DW_V':'TP_DW_V',
#'TP_DW_H':'TP_DW_H','TOI_FLSH':'TOI_flsh','TOI_PUMP':'TOI_pump',
#'TOI_READ':'TOI_read','TOI_CINJ':'TOI_CInj',
#'INVFLSHP':'Invflshp','INVFLUSH':'Invflush','FLUSHES':'Flushes',
#'VSTART':'Vstart','VEND':'Vend',
#'OVRSCN_H':'Ovrscn_H','CLK_ROE':'CLK_ROE','CNVSTART':'CnvStart',
#'SUMWELL':'SumWell','INISWEEP':'IniSweep','SPW_CLK':'SPW_clk',
#'FPGA_VER':'FPGA_ver','EGSE_VER':'EGSE_ver',
#'M_STEPS':'M_Steps','M_ST_SZE':'M_st_Sze',
#'WAVELENG':'Wavelength','MIRR_POS':'Mirr_pos',
#'CHMB_PRE':'Chmb_pre','CCD1_SN':'CCD1_SN','CCD2_SN':'CCD2_SN','CCD3_SN':'CCD3_SN',
#'ROE_SN':'ROE_SN','CALSCRPT':'CalScrpt','COMMENTS':None,
#'TMPCCD1T':'R1CCD1TT','TMPCCD1B':'R1CCD1TB',
#'TMPCCD2T':'R1CCD2TT','TMPCCD2B':'R1CCD2TB',
#'TMPCCD3T':'R1CCD3TT','TMPCCD3B':'R1CCD3TB',
#'IDL_V':'IDL_V','IDH_V':'IDH_V',
#'IG1_T1_V':'IG1_T1_V','IG1_T2_V':'IG1_T2_V','IG1_T3_V':'IG1_T3_V',
#'IG1_B1_V':'IG1_B1_V','IG1_B2_V':'IG1_B2_V','IG1_B3_V':'IG1_B3_V',
#'IG2_T_V':'IG2_T_V','IG2_B_V':'IG2_B_V',
#'OD_T1_V':'OD_T1_V','OD_T2_V':'OD_T2_V','OD_T3_V':'OD_T3_V',
#'OD_B1_V':'OD_B1_V','OD_B2_V':'OD_B2_V','OD_B3_V':'OD_V3_V',
#'RD_T_V':'RD_T_V','RD_B_V':'RD_B_V'}}


class Extension():
    """Extension Class"""
    
    def __init__(self,data,header=None,label=None,headerdict=None):
        """ """
        
        self.data = data
        
        if header is None:
            header = fts.Header()
        
        if headerdict is not None:
            for key in headerdict:
                header[key] = headerdict[key]
                
        self.header = header
        
        if label is not None:
            self.header['EXTNAME'] = label
        
        self.label = label



class CCD(object):
    """Class of CCD objects. 
    Euclid Images as acquired by ELVIS software (Euclid LabView Imaging Software).
    
    
    The class has been extended to handle multi-extension images. This is useful
    to also "host" calibration data-products, such as Flat-Fields.
    
    """
    
    def __init__(self,infits=None,extensions=[-1],getallextensions=False):
        """ """
        
        self.extnames = []
        self.extensions = []
        
        
        if infits is not None:
        
            assert type(infits) is str, "'%s' can't be a name for a file!" % infits
            assert isthere(infits), 'infits is just not there :-('
                        
            self.loadfromFITS(infits,extensions,getallextensions)

        else:
            
            self.extensions = []
            self.extnames = []
        
        self.nextensions = len(self.extensions)
        
        self.NAXIS1 = NAXIS1
        self.NAXIS2 = NAXIS2        
        self.shape = (NAXIS1,NAXIS2)
        
        for iext in range(self.nextensions):
            if self.extensions[iext].data is not None:
                assert self.shape == self.extensions[iext].data.shape                
        
        self.prescan = prescan
        self.overscan = overscan
        self.gain = dict(E=3.1,F=3.1,G=3.1,H=3.1)
        self.rn = dict(E=4.5,F=4.5,G=4.5,H=4.5)
        
        self.QuadBound = QuadBound 
        
        self.masked = False
    
    
    def loadfromFITS(self,fitsfile,extensions=[-1],getallextensions=False):
        
        hdulist = fts.open(fitsfile)
        
        
        nextensions = len(hdulist)
            
        if getallextensions:
            extensions = np.arange(nextensions)
            
        for iext in extensions:
                
            hdu = hdulist[iext]
            
            
            if hdu.data is not None:
                data = hdu.data.transpose().astype('float32').copy()
            else: data = None
            header = hdu.header
                
            if 'EXTNAME' in hdu.header:
                label = hdu.header['EXTNAME']
            else:
                label = None
                
            self.extensions.append(Extension(data,header,label))
            self.extnames.append(label)
            
        hdulist.close()

    
    def add_extension(self,data,header=None,label=None,headerdict=None):
        """ """    
        self.extensions.append(Extension(data,header,label,headerdict))
        self.nextensions += 1
        
    
    def get_quad(self,Quadrant,canonical=False,extension=-1):
        """Returns a quadrant in canonical or non-canonical orientation.
        
        :param Quadrant: Quadrant, one of 'E', 'F', 'G', 'H'
        :type Quadrant: char
        
        :param canonical: 
        
        Canonical [True] = with readout-node at pixel index (0,0) regardless of quadrant. This is the orientation which corresponds to the data-reading order (useful for cross-talk measurements, for example).
        Non-Canonical [False] = with readout-node at corner matching placement of quadrant on the CCD. This is the orientation that would match the representation of the image on DS9.        

        :type canonical: bool
        
        :param extension: extension number. Default = -1 (last)
        :type extension: int
        
        """
        
        edges = self.QuadBound[Quadrant]        
        Qdata = self.extensions[extension].data[edges[0]:edges[1],edges[2]:edges[3]]
        
        if canonical:
            if Quadrant == 'E': Qdata = Qdata[:,::-1].copy()
            elif Quadrant == 'F': Qdata = Qdata[::-1,::-1].copy()
            elif Quadrant == 'G': Qdata = Qdata[::-1,:].copy()
            elif Quadrant == 'H': Qdata = Qdata[:,:].copy()
        
        
        return Qdata
        
    def get_cutout(self,corners,Quadrant,canonical=False,extension=-1):
        """Returns a cutout from the CCD image, either in 
        canonical or non-canonical orientation.
        
        
        :param corners: [x0,x1,y0,y1]
        :type corners: list (of int)
        
        :param Quadrant: Quadrant, one of 'E', 'F', 'G', 'H'
        :type Quadrant: char
        
        :param canonical: 
         Canonical [True] = with readout-node at pixel index (0,0) regardless of quadrant. 
         This is the orientation which corresponds to the data-readin order 
         (useful for cross-talk measurements, for example).
         Non-Canonical [False] = with readout-node at corner matching placement of  
         quadrant on the CCD. This is the orientation that would match the representation 
         of the image on DS9.        
        :type canonical: bool
        
        :param extension: extension number. Default = -1 (last)
        :type extension: int
        
        
        """       
        Qdata = self.get_quad(Quadrant,canonical,extension)
        cutout = Qdata[corners[0]:corners[1],corners[2]:corners[3]].copy()
        return cutout
        
        
    def set_quad(self,inQdata,Quadrant,canonical=False,extension=-1):
        """ """
        edges = self.QuadBound[Quadrant]
        #Qdata = self.data[edges[0]:edges[1],edges[2]:edges[3]]
        
        if canonical:
            if Quadrant == 'E': Qdata = inQdata[:,::-1].copy()
            elif Quadrant == 'F': Qdata = inQdata[::-1,::-1].copy()
            elif Quadrant == 'G': Qdata = inQdata[::-1,:].copy()
            elif Quadrant == 'H': Qdata = inQdata[:,:].copy()
        else:
            Qdata = inQdata.copy()
        
        self.extensions[extension].data[edges[0]:edges[1],edges[2]:edges[3]] = Qdata.copy()
        
        return None


    def getsectioncollims(self,QUAD):
        """Returns limits of sections: prescan, image and overscan"""
        
        semiNAXIS1 = self.NAXIS1/2
        
        if QUAD in ['E','H']:
                
            prestart = 0
            preend = self.prescan-1
            ovstart = semiNAXIS1-self.overscan 
            ovend = ovstart + self.overscan - 1
            imgstart = preend+1
            imgend = ovstart-1
            
        elif QUAD in ['F','G']:
            
            ovstart = 0
            ovend = self.overscan-1
            prestart = semiNAXIS1-self.prescan
            preend = prestart + self.prescan-1
            imgstart = ovend+1
            imgend = prestart-1
                
        return (prestart,preend,imgstart,imgend,ovstart,ovend)


   
    def get_stats(self,Quadrant,sector='img',statkeys=['mean'],trimscan=[0,0],
                  ignore_pover=True,extension=-1):
        """ """
        
        Qdata = self.get_quad(Quadrant,canonical=True,extension=extension)
        
        if isinstance(Qdata,np.ma.masked_array):
            stat_dict = dict(mean=np.ma.mean,median=np.ma.median,std=np.ma.std)
        else:
            stat_dict = dict(mean=np.mean,median=np.median,std=np.std)
        
        if ignore_pover:
            vlims = [0,NrowsCCD]
        else:
            vlims = [0,None]
        
        if sector == 'pre':
            hlims = [0,self.prescan]
        elif sector == 'ove':
            hlims = [self.NAXIS1/2-self.overscan,self.NAXIS1/2]
        elif sector == 'img':
            hlims = [self.prescan,self.NAXIS1/2-self.overscan]
        elif sector == 'all':
            hlims = [0,None]
            
        hlims[0] += trimscan[0]
        hlims[1] -= trimscan[1]
        
        results = []
        
        for statkey in statkeys:
            
            results.append(stat_dict[statkey](Qdata[hlims[0]:hlims[1],vlims[0]:vlims[1]]))
        
        return results
    
    def sub_offset(self,Quad,method='row',scan='pre',trimscan=[3,2],
                   ignore_pover=True,extension=-1):
        """ """
        
        if self.masked:
            median= np.ma.median
        else:
            median = np.median
        
        quaddata = self.get_quad(Quad,canonical=True,extension=extension)
        
        
        if scan == 'pre':
            hlims = [0,self.prescan]
        elif scan == 'ove':
            hlims = [self.NAXIS1/2-self.overscan,self.NAXIS1/2]
        else:
            sys.exit('ccd.sub_offset: scan=%s unkonwn' % scan)
            
        hlims[0] += trimscan[0]
        hlims[1] -= trimscan[1]
        
        if ignore_pover:
            vlims = [0,NrowsCCD]
        else:
            vlims = [0,None]

        if method == 'row':
            offsets = []
            for ix in range(self.NAXIS2/2):
                offset = median(quaddata[hlims[0]:hlims[1],ix])
                if self.masked : offset = offset.data
                quaddata[:,ix] -= offset
                offsets.append(offset)
            
        elif method == 'median':
            
            offset = median(quaddata[hlims[0]:hlims[1],vlims[0]:vlims[1]])
            #if self.masked : offset = offset.data
            quaddata -= offset
            offsets = [offset]
        
        B = self.QuadBound[Quad]
        self.extensions[extension].data[B[0]:B[1],B[2]:B[3]] = self.flip_tocanonical(quaddata,Quad).copy()
        

        return offsets
        
    
    def sub_bias(self,superbias,extension=-1):
        """Subtracts a superbias"""
        
        assert self.shape == superbias.shape
        self.extensions[extension].data -= superbias
        
    
    def divide_by_flatfield(self,FF,extension=-1):
        """Divides by a Flat-field"""
        print 'TODO: ccd.CCD.divide_by_flatfield needs improvements: handling of masked values'
        assert self.shape == FF.shape
        self.extensions[extension].data /= FF
        
    
    def flip_tocanonical(self,array,Quad):
        
        if Quad == 'E': return array[:,::-1].copy()
        elif Quad == 'F': return array[::-1,::-1].copy()
        elif Quad == 'G': return array[::-1,:].copy()
        elif Quad == 'H': return array.copy()
    
    def do_Vscan_Mask(self,VSTART,VEND):
        """ """
        
        VscanMask = np.ones((self.NAXIS1,self.NAXIS2),dtype='bool')
        
        for Quad in self.QuadBound.keys():
            
            B = self.QuadBound[Quad]
            
            tmp = self.flip_tocanonical(VscanMask[B[0]:B[1],B[2]:B[3]],Quad)
            tmp[:,VSTART:VEND+1] = False
            VscanMask[B[0]:B[1],B[2]:B[3]] = self.flip_tocanonical(tmp,Quad).copy()
        
        
        return VscanMask
    
    def get_mask(self,mask):
        """ """
        assert self.shape == mask.shape
        
        for iext in range(self.nextensions):
            
            masked = np.ma.masked_array(self.extensions[iext].data,mask)
            self.extensions[iext].data = masked.copy()
        
        self.masked = True
        
    
    def writeto(self,fitsf,clobber=False,unsigned16bit=False):
        """ """
        
        #prihdu = fts.PrimaryHDU()
        
        firstextension = self.extensions[0]
          
        if firstextension.data is not None:
            if self.masked: pridata = firstextension.data.data.transpose().copy()
            else: pridata = firstextension.data.transpose().copy()
        else:
            pridata = None
            
        prihdr = firstextension.header
                
        prihdu = fts.PrimaryHDU(data=pridata,header=prihdr)

          
        comments = ['FITS file generated by vison.datamodel.ccd',
        'point of contact: Ruyman Azzollini (r.azzollini_at_ucl.ac.uk)',]
        for comm in comments:
            prihdu.header.add_comment(comm)
            
        hdulist = fts.HDUList([prihdu])
        
        if self.nextensions > 1:
            
            for iext in range(1,self.nextensions):
                
                if self.masked: idata = self.extensions[iext].data.data.transpose().copy()
                else: idata = self.extensions[iext].data.transpose().copy()
                
                iheader = self.extensions[iext].header
                iname = self.extensions[iext].label
                
                ihdu = fts.ImageHDU(data=idata,header=iheader,name=iname)
                
                if unsigned16bit:
                    ihdu.scale('int16', '', bzero=32768)
                    ihdu.header.add_history('Scaled to unsigned 16bit integer!')
            
                hdulist.append(ihdu)
        
        hdulist.writeto(fitsf,overwrite=clobber)
    
    
    def simadd_flatilum(self,levels=dict(E=0.,F=0.,G=0.,H=0.),extension=-1):
        """ """

        for Q in Quads:
            quaddata = self.get_quad(Q,canonical=True,extension=extension)
            quaddata[self.prescan:-self.overscan,0:NrowsCCD] += levels[Q]
            self.set_quad(quaddata,Q,canonical=True,extension=extension)
    
    
    def simadd_points(self,flux,fwhm,CCDID='CCD1',dx=0,dy=0,extension=-1):
        """ """
        from vison.point.lib import Point_CooNom
        from vison.point.models import fgauss2D
        
        sigma = fwhm/2.355
        i0 = flux/(2.*np.pi*sigma**2.)
        
        nx = 15
        ny = 15

        
        for Q in Quads:
            quaddata = self.get_quad(Q,canonical=False,extension=extension).copy()
            
            point_keys = Point_CooNom[CCDID][Q].keys()
            
            B = self.QuadBound[Q]
            
            for pkey in point_keys:
                xp,yp = Point_CooNom[CCDID][Q][pkey]
                
                x0 = xp - B[0] + dx
                y0 = yp - B[2] + dy
                
                xmin = int(np.round(x0))-nx/2
                xmax = xmin+nx
                ymin = int(np.round(y0))-ny/2
                ymax = ymin+ny
                
                x = np.arange(xmin,xmax)
                y = np.arange(ymin,ymax)
                
                xx,yy = np.meshgrid(x,y,indexing='xy')
                
                p = [i0,x0,y0,sigma,sigma,0.]
                point_stamp = fgauss2D(xx,yy,p)
                
                quaddata[xmin:xmax,ymin:ymax] += point_stamp.copy()
                
            self.set_quad(quaddata,Q,canonical=False,extension=extension)
        
        
        
    
    def simadd_bias(self,levels=dict(E=2000,F=2000,G=2000,H=2000),extension=-1):
        
         for Q in Quads:            
            B = self.QuadBound[Q]
            self.extensions[extension].data[B[0]:B[1],B[2]:B[3]] += levels[Q]
    
    
    def simadd_ron(self,extension=-1):
        """ """
        
        noise = np.random.normal(loc=0.0, scale=self.rn['E'],size=self.extensions[extension].data.shape)
        self.extensions[extension].data += noise
        
    
    def simadd_poisson(self,extension=-1):
        """ """
        
        gain = self.gain['E']
        
        image_e = self.extensions[extension].data.copy()*gain
        
        rounded = np.rint(image_e)
        residual = image_e - rounded #ugly workaround for multiple rounding operations...
        rounded[rounded < 0.0] = 0.0
        image_e = np.random.poisson(rounded).astype(np.float64)
        image_e += residual
        
        image = image_e / gain
        
        self.extensions[extension].data = image.copy()
            
    def sim_window(self,vstart,vend,extension=-1):
        """ """
        
        mask = np.ones((self.NAXIS1/2,self.NAXIS2/2),dtype='int8')
        mask[:,vstart:vend] = 0
        
        for Q in Quads:            
            qdata = self.get_quad(Q,canonical=True,extension=extension)
            qdata[np.where(mask)] = 0
            self.set_quad(qdata,Q,canonical=True,extension=extension)
        
    def simadd_injection(self,levels,on=NrowsCCD,off=0,extension=-1):
        
        self.simadd_flatilum(levels=levels,extension=-1)
        
        # ON/OFF masking
        
        NX,NY = self.NAXIS1/2,self.NAXIS2/2
        
        mask_onoff = np.zeros((NX,NY),dtype='int8')
        
        y0 = 0
        Fahrtig = False
        
        while not Fahrtig:
            y1 = min(y0 + on,NY)
            mask_onoff[:,y0:y1] = 1
            y0 = y1 + off
            if (y0>NY-1):
                break

        #from astropy.io import fits as fts
        #fts.writeto('mask.fits',mask_onoff.transpose().astype('int32'),overwrite=True)
        
        for Q in Quads:
            qdata = self.get_quad(Q,canonical=True,extension=extension)
            qdata[np.where(mask_onoff==0)] = 0.
            self.set_quad(qdata,Q,canonical=True,extension=extension)
        
    
def test_create_from_scratch():
    """ """
    
    NAXIS1,NAXIS2 = 4238,4172
    
    img = np.ones((NAXIS1,NAXIS2),dtype='float32')
    eimg = np.ones((NAXIS1,NAXIS2),dtype='float32') * 0.1
    
    ccdout = CCD()
    
    fitsname = 'test_create_from_scratch.fits'
    
    ccdout.add_extension(data=img,label='IMAGE')
    ccdout.add_extension(data=eimg,label='UNCERTAINTY')
    ccdout.writeto(fitsname,clobber=True)
    
    ccdin = CCD(fitsname,getallextensions=True)
    
    print 'Number of extensions = %i' % ccdin.nextensions
    stop()

def test_load_ELVIS_fits():
    """ """
    
    fitsname = '/home/raf/WORK/EUCLID/REPOS/vison/vison/data/EUC_2112_231016D_135042T_ROE1_CCD1.fits'
    
    ccd = CCD(infits=fitsname,getallextensions=True)
    
    ccd.writeto('ccd_test_load_ELVIS_fits.fits',clobber=True)
    
    
if __name__ == '__main__':
    
    #test_create_from_scratch()
    test_load_ELVIS_fits()
    