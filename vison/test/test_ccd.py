#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Unit-testing for CCD class.


Created on Mon May  7 09:47:07 2018

:author: Ruyman Azzollini

"""
# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import unittest
from vison.datamodel import ccd
# END IMPORT

class TestCCDClass(unittest.TestCase):
    
    def setUp(self):
        self.NAXIS1 = 4238
        self.wQ = self.NAXIS1/2
        self.NAXIS2 = 4172
        self.hQ = self.NAXIS2/2
        self.offset = 1000.
        self.badpixel = (100,100) # bad pixel on each Quad, canonical orientation
        self.prescan = ccd.prescan
        self.overscan = ccd.overscan
        self.tolerance = 1.e-7
        self.gain = 3.5
        self.RON = 4.5 / self.gain
        
        img = np.zeros((self.NAXIS1,self.NAXIS2),dtype='float32') 
        eimg = np.ones((self.NAXIS1,self.NAXIS2),dtype='float32') * self.RON
                
        self.ccdobj = ccd.CCD()
        self.ccdobj.add_extension(data=img,label='IMAGE')
        self.ccdobj.add_extension(data=eimg,label='UNCERTAINTY')
        
        NX = self.NAXIS1/2
        NY = self.NAXIS2/2
        
        for Q in self.ccdobj.Quads:
        
            
            Y, _ = np.meshgrid(np.arange(NY), np.arange(NX), indexing='xy')
        
            Qramp = Y.copy()*1.
            Qramp[0:self.prescan,:] = 0.
            Qramp[-self.overscan:,:] = 0.
            Qramp += self.offset
            
            Qramp[self.badpixel] = np.nan
        
            self.ccdobj.set_quad(Qramp.copy(),Q,canonical=True,extension=0)
        
        mask = np.isnan(self.ccdobj.extensions[0].data)
        self.ccdobj.get_mask(mask)
        
    
    def test_ccdobj_has_extensions(self):
        self.assertEqual(self.ccdobj.nextensions,2,
            msg='ccdobj has %i extensions (expected 2)' % self.ccdobj.nextensions)
        
    def test_ccdobj_has_right_shape(self):
        self.assertEqual(self.ccdobj.shape,(self.NAXIS1,self.NAXIS2),
                         msg='shape=%s, expected=%s' %\
                         (self.ccdobj.shape.__repr__(),(self.NAXIS1,self.NAXIS2).__repr__()))
    
    def test_get_mask(self):
        """ """
        mask = np.isnan(self.ccdobj.extensions[0].data)
        self.ccdobj.get_mask(mask)
        self.assertTrue(self.ccdobj.masked)
    
    def _validate_coo_conversion(self,f,incoos,outcoos):
        
        for i in range(len(incoos)):
            point = incoos[i]
            xin,yin,Q = point
            xou,you = f(xin,yin,Q)
            dist = ((outcoos[i][0]-xou)**2.+(outcoos[i][1]-you)**2.)**0.5
            self.assertAlmostEqual(dist,0.,delta=self.tolerance,
                             msg='running subtest %i' %  i)
    
    def test_cooconv_Qrel_2_CCD(self):
        
        incoos = [(100,200,'E'),
                  (100,200,'F'),
                  (100,200,'G'),
                  (100,200,'H')]
        outcoos = [(incoos[0][0],incoos[0][1]+self.hQ),
                   (incoos[1][0]+self.wQ,incoos[1][1]+self.hQ),
                   (incoos[2][0]+self.wQ,incoos[2][1]),
                   (incoos[3][0],incoos[3][1])]
        
        self._validate_coo_conversion(self.ccdobj.cooconv_Qrel_2_CCD,
                                      incoos,outcoos)
    
    def test_cooconv_Qcan_2_CCD(self):
        
        hQ = self.NAXIS2/2
        wQ = self.NAXIS1/2
        
        incoos = [(100,100,'E'),
                  (100,100,'F'),
                  (100,100,'G'),
                  (100,100,'H')]
        outcoos = [(incoos[0][0],2*hQ-incoos[0][1]-1.),
                   (2.*wQ-incoos[1][0]-1.,2*hQ-incoos[1][1]-1.),
                   (2.*wQ-incoos[2][0]-1.,incoos[2][1]),
                   (incoos[3][0],incoos[3][1])]

        self._validate_coo_conversion(self.ccdobj.cooconv_Qcan_2_CCD,
                                      incoos,outcoos)
        
    def test_cooconv_CCD_2_Qrel(self):
        
        hQ = self.NAXIS2/2
        wQ = self.NAXIS1/2
        
        incoos = [(100,100+hQ,'E'),
                  (100+wQ,100+hQ,'F'),
                  (100+wQ,100,'G'),
                  (100,100,'H')]
        
        outcoos = [(100.,100.),
                   (100.,100.),
                   (100.,100.),
                   (100.,100.)]

        self._validate_coo_conversion(self.ccdobj.cooconv_CCD_2_Qrel,
                                      incoos,outcoos)

        
    def test_cooconv_CCD_2_Qcan(self):
        
        hQ = self.NAXIS2/2
        wQ = self.NAXIS1/2
        
        incoos = [(100,2*hQ-100-1.,'E'),
                  (2*wQ-100.-1.,2.*hQ-100-1.,'F'),
                  (2.*wQ-100-1.,100,'G'),
                  (100,100,'H')]
        
        outcoos = [(100.,100.),
                   (100.,100.),
                   (100.,100.),
                   (100.,100.)]

        self._validate_coo_conversion(self.ccdobj.cooconv_CCD_2_Qcan,
                                      incoos,outcoos)

        
    def test_cooconv_Qrel_2_Qcan(self):
        
        hQ = self.NAXIS2/2
        wQ = self.NAXIS1/2
        
        incoos = [(100,hQ-100-1.,'E'),
                  (wQ-100.-1.,hQ-100-1.,'F'),
                  (wQ-100-1.,100,'G'),
                  (100,100,'H')]
        
        outcoos = [(100.,100.),
                   (100.,100.),
                   (100.,100.),
                   (100.,100.)]

        self._validate_coo_conversion(self.ccdobj.cooconv_Qrel_2_Qcan,
                                      incoos,outcoos)

            
    def test_cooconv_Qcan_2_Qrel(self):
        hQ = self.NAXIS2/2
        wQ = self.NAXIS1/2
        
        incoos = [(100.,100.,'E'),
                   (100.,100.,'F'),
                   (100.,100.,'G'),
                   (100.,100.,'H')]
        
        outcoos = [(100,hQ-100-1.),
                  (wQ-100.-1.,hQ-100-1.),
                  (wQ-100-1.,100),
                  (100,100)]

        self._validate_coo_conversion(self.ccdobj.cooconv_Qrel_2_Qcan,
                                      incoos,outcoos)
        
    
    #def test_dummyrebin(self):
    #    pass
    
    def test_get_tile_coos(self):
        
        wpx = 200
        hpx = 150
        Qs = self.ccdobj.Quads
        
        def get_img_bounds(Q):
            _, _, _imgstarth, _imgendh, _, _ = self.ccdobj.getsectioncollims(Q)
            _imgstartv, _imgendv, _, _ = self.ccdobj.getsectionrowlims(Q)
            return (_imgstarth,_imgendh, _imgstartv, _imgendv)
        
        def get_urpix(llpix,wpx,hpx):
            Nsamps = len(llpix)
            urpix = [(llpix[i][0]+wpx-1,llpix[i][1]+hpx-1) for i in range(Nsamps)]
            return urpix
        
        def are_within_image(llpix,urpix, wpx, hpx, img_bounds):
            arewithinimage = [(llpix[i][0]>=img_bounds[0]) & (urpix[i][0] <= img_bounds[1]) &\
                             (llpix[i][1]>=img_bounds[2]) & (urpix[i][1] <= img_bounds[3]) \
                             for i in range(len(llpix))]
            return arewithinimage
        
        for i,Q in enumerate(Qs):
            img_bounds = get_img_bounds(Q)           
            itiles = self.ccdobj.get_tile_coos(Q,wpx,hpx)
            ccpix = itiles['ccpix']
            llpix = itiles['llpix']
            Nsamps = itiles['Nsamps']
            
            comm_msg = 'Failed on checks on Q=%s' % Q
            
            self.assertTrue((wpx == itiles['wpx']) and (hpx == itiles['hpx']),msg='check 1, %s' % comm_msg)
            self.assertTrue(len(ccpix) == len(llpix) == Nsamps, msg='check 2, %s' % comm_msg)
            
            urpix = get_urpix(llpix,wpx,hpx)
            arewithinimage = are_within_image(llpix, urpix, wpx, hpx, img_bounds)
            
            self.assertTrue(np.all(arewithinimage), msg='check 3, %s' % comm_msg)
            
            areas = [(urpix[i][0]-llpix[i][0]+1)*(urpix[i][1]-llpix[i][1]+1) \
                     for i in range(Nsamps)]
            
            
            self.assertTrue(np.all(np.isclose(np.array(areas), wpx*hpx)), msg='check 4, %s' % comm_msg)
            
    
    def test_get_tiles_stats(self):
        
        Q = 'E'
        wpx = 250
        hpx = 150
        tile_coos = self.ccdobj.get_tile_coos(Q,wpx,hpx)
        statkey = 'mean'
        
        res = self.ccdobj.get_tiles_stats(Q, tile_coos, statkey, extension=-1)
        
        self.assertTrue(len(res) == tile_coos['Nsamps'])
        self.assertTrue(np.all(np.isclose(res,self.RON)))
        
    
    def test_get_cutout(self):
        prescan = self.ccdobj.prescan
        corners = [prescan,prescan+2,0,2]
        expected = np.array([[self.offset,self.offset+1.],
                             [self.offset,self.offset+1.]])
        for iQ,Q in enumerate(self.ccdobj.Quads):
            _cut = self.ccdobj.get_cutout(corners,Q,canonical=True,extension=0)
            
            self.assertAlmostEqual(np.abs(_cut-expected).sum(),0.,delta=self.tolerance,
                                   msg='Failed at check %i,Q=%s' % (iQ,Q))
    
    def test_getsectioncollims(self):
        pre = self.ccdobj.prescan
        over = self.ccdobj.overscan
        img = self.wQ-pre-over
        expected = dict(E=(0, pre-1, pre, img+pre-1, 
                           img+pre, img+pre+over-1),
                        F=(img+over, img+over+pre-1, 
                           over, img+over-1, 0, over-1),
                        G=(img+over, img+over+pre-1, 
                           over, img+over-1, 0, over-1),
                        H=(0, pre-1, pre, 
                           img+pre-1, img+pre, img+pre+over-1))
        for iQ,Q in enumerate(self.ccdobj.Quads):
            res = np.array(self.ccdobj.getsectioncollims(Q))
            shouldbe0 = np.abs(np.array(expected[Q]) - res).sum()
            sizesareok = (((res[1]-res[0]+1)==pre) &\
                          ((res[3]-res[2]+1)==img) &\
                          ((res[5]-res[4]+1)==over))
            self.assertTrue(np.isclose(shouldbe0,0.,rtol=self.tolerance), 
                            msg='values error... Failed at check %i,Q=%s' % (iQ,Q))
            self.assertTrue(sizesareok,msg='dimensions error... Failed at check %i,Q=%s' % (iQ,Q))
            
    
    def test_getsectionrowlims(self):
        
        vover = self.ccdobj.voverscan
        img = self.ccdobj.NrowsCCD
        expected = dict(E=(vover, vover+img-1, 0, vover-1),
                        F=(vover, vover+img-1, 0, vover-1),
                        G=(0, img-1, img, img+vover-1),
                        H=(0, img-1, img, img+vover-1))
        
        for iQ,Q in enumerate(self.ccdobj.Quads):
            res = np.array(self.ccdobj.getsectionrowlims(Q))
            shouldbe0 = np.abs(np.array(expected[Q]) - res).sum()
            sizesareok = ((res[1]-res[0]+1)==img) & \
                         ((res[3]-res[2]+1)==vover)
                         
            self.assertTrue(np.isclose(shouldbe0,0.,rtol=self.tolerance), 
                            msg='values error... Failed at check %i,Q=%s' % (iQ,Q))
            self.assertTrue(sizesareok,msg='dimensions error... Failed at check %i,Q=%s' % (iQ,Q))
            
    
    def test_get_stats(self):
        Nimg = self.ccdobj.NrowsCCD
        offset = self.offset
        pre=self.ccdobj.prescan
        over= self.ccdobj.overscan
        NY = self.NAXIS2/2
        NX = self.NAXIS1/2
        Y, _ = np.meshgrid(np.arange(NY,dtype='float32'), 
                           np.arange(NX,dtype='float32'), indexing='xy')
        Y = Y.astype('float32') + offset
        Y[self.badpixel] = np.nan
        mask = np.isnan(Y)
        Y = np.ma.masked_array(Y,mask)
        expected = np.ma.mean(Y[pre:-over,0:Nimg])
        
        for iQ,Q in enumerate(self.ccdobj.Quads):
            res = self.ccdobj.get_stats(Q,sector='img',statkeys=['mean'],
                                        trimscan=[0,0],ignore_pover=True,
                                        extension=0)
            self.assertAlmostEqual(res[0],expected,delta=self.tolerance,
                                   msg='Failed at check %i,Q=%s' % (iQ,Q))
        
    
    def test_sub_offset(self):
        
        expected = self.offset
        
        for iQ,Q in enumerate(self.ccdobj.Quads):
            
            res = self.ccdobj.sub_offset(Q,method='median',scan='pre',
                                        trimscan=[0,0],ignore_pover=True,
                                        extension=0)
            self.assertAlmostEqual(res[0],expected,delta=self.tolerance,
                                   msg='Failed at check %i,Q=%s' % (iQ,Q))
        
    
    def test_sub_bias(self):
        superbias = np.zeros_like(self.ccdobj.extensions[0].data)+self.offset
        self.ccdobj.sub_bias(superbias,extension=0)
        res = self.ccdobj.get_stats('E',sector='pre',statkeys=['mean'],
                                        trimscan=[0,0],ignore_pover=True,
                                        extension=0)
        self.assertAlmostEqual(res[0],0.,delta=self.tolerance)
        
    
    def test_divide_by_flatfield(self):
        FF = np.ones_like(self.ccdobj.extensions[0].data)*2.
        img_bef = self.ccdobj.extensions[0].data.copy()
        self.ccdobj.divide_by_flatfield(FF,extension=0)
        img_aft = self.ccdobj.extensions[0].data.copy()
        self.assertAlmostEqual(np.ma.mean((img_aft/img_bef)),0.5,
                               delta=self.tolerance)
    
    def test_get_1Dprofile(self):
        
        expected = np.arange(0,self.ccdobj.NrowsCCD,dtype='float32')+self.offset
        
        for iQ,Q in enumerate(self.ccdobj.Quads):
                            
            res = self.ccdobj.get_1Dprofile('E',orient='ver',area='img',
                                            stacker='median',vstart=0,
                                            vend=self.ccdobj.NrowsCCD,
                                            extension=0)
            resy = res.data['y']
            
            self.assertAlmostEqual(np.abs(expected-resy).sum(),0.,
                                   delta=self.tolerance,
                                   msg='Failed at check %i,Q=%s' % (iQ,Q))
        
    
    def test_get_region2Dmodel(self):
        
        for iQ,Q in enumerate(self.ccdobj.Quads):
                            
            res = self.ccdobj.get_region2Dmodel('E',area='img',kind='poly2D',
                        pdegree=1,doFilter=False,doBin=False,
                        vend=self.ccdobj.NrowsCCD,canonical=True,
                        extension=0)
            rms = np.sqrt((res.imgmodel-res.img)**2.).mean()
            
            self.assertAlmostEqual(rms,0.,
                                   delta=self.tolerance,
                                   msg='Failed at check %i,Q=%s' % (iQ,Q))            
        
    
    def test_extract_region(self):
        
        Quads = self.ccdobj.Quads
        
        subreg_r, BB_r = self.ccdobj.extract_region(Quads[0],area='img',vstart=0,vend=2086,
                                                Full=True,canonical=True,extension=0)    
        
        for Q in Quads[1:]:
        
            subreg, BB = self.ccdobj.extract_region(Q,area='img',vstart=0,vend=2086,
                                   Full=True,canonical=True,extension=0)
            self.assertAlmostEqual(np.abs(subreg-subreg_r).sum(),0.,delta=self.tolerance,
                                   msg='Failed at Q=%s' % Q)
    

 
if __name__ == '__main__':
    unittest.main()