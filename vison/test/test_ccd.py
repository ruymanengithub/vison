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
        
        img = np.zeros((self.NAXIS1,self.NAXIS2),dtype='float32') 
        eimg = np.ones((self.NAXIS1,self.NAXIS2),dtype='float32') * 1.4
                
        self.ccdobj = ccd.CCD()
        self.ccdobj.add_extension(data=img,label='IMAGE')
        self.ccdobj.add_extension(data=eimg,label='UNCERTAINTY')
        
        NX = self.NAXIS1/2
        NY = self.NAXIS2/2
        Y, _ = np.meshgrid(np.arange(NY), np.arange(NX), indexing='xy')
        
        Qramp = Y.copy()*1.
        Qramp[0:self.prescan,:] = 0.
        Qramp[-self.overscan:,:] = 0.
        Qramp += self.offset
        
        Qramp[self.badpixel] = np.nan
        
        for Q in self.ccdobj.Quads:
            self.ccdobj.set_quad(Qramp.copy(),Q,canonical=False,extension=0)
        
    
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
        pass
    
    def test_get_tiles_stats(self):
        pass
    
    def test_get_cutout(self):
        pass
    
    def test_getsectioncollims(self):
        pass
    
    def test_getsectionrowlims(self):
        pass
    
    def test_get_stats(self):
        pass
    
    def test_sub_offset(self):
        pass
    
    def test_sub_bias(self):
        pass
    
    def test_divide_by_flatfield(self):
        pass
    
    def test_get_1Dprofile(self):
        pass
    
    def test_get_region2Dmodel(self):
        pass
    
    def test_extract_region(self):
        pass
    

 
if __name__ == '__main__':
    unittest.main()