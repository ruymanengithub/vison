#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Sextractor interface.

Created on Thu May 17 13:29:05 2018

@author: raf
"""

# IMPORT STUFF
from pdb import set_trace as stop
import astromatic_wrapper as aw
import os
import string as st

from vison.datamodel import ccd as ccdmod
# END IMPORT

default_params = ['NUMBER','EXT_NUMBER','X_IMAGE','Y_IMAGE','A_IMAGE','B_IMAGE','ELONGATION','FWHM_IMAGE','MAG_AUTO']

class VSExtractor(object):
    
    
    def __init__(self,FITS=None,ccdobj=None):
        """ """        
        self.files = dict()
        
        self.internals = dict(
                MINAREA=5,
                MAG_ZEROPOINT=20.,
                SATUR_LEVEL=65535,
                SEEING_FWHM = 1.5,
                PIXEL_SCALE = 1.,
                GAIN = 1.
                )
        self.internals['params'] = default_params
    
    def update_internals(self,inputs):
        assert isinstance(inputs,dict)
        self.internals.update(inputs)
    
    def get_sex_kwargs(self,catroot,config=None,checks=None):
        
        catname = '%s.ldac' % catroot
        
        configdefaults = dict(
                CATALOG_NAME = catname,
                DETECT_MINAREA='%i' % self.internals['MINAREA'],
                DETECT_THRESH='15.,%.1f' % self.internals['MAG_ZEROPOINT'],
                ANALYSIS_THRESH='15.,%.1f' % self.internals['MAG_ZEROPOINT'],
                FILTER='N',
                SATUR_LEVEL='%i' % self.internals['SATUR_LEVEL'],
                MAG_ZEROPOINT='%.1f' % self.internals['MAG_ZEROPOINT'],
                GAIN='%.1f' % self.internals['GAIN'],
                PIXEL_SCALE='%.1f' % self.internals['PIXEL_SCALE'],
                SEEING_FWHM='%.1f' % self.internals['SEEEING_FWHM'],
                BACKPHOTO_TYPE='LOCAL',
                BACK_SIZE='500',
                BACK_FILTERSIZE='3')

        if checks is not None:
            assert isinstance(checks,list)
            assert isinstance(checks[0],str)
            
            checks = [item.upper() for item in checks]
            
            checkimage_name = st.join(['%s_%s.fits' % \
                                       (catroot,check[0:4]) for check in checks],',')
            
            checkimage_type = st.join(checks,',')
            
            configdefaults.update(dict(
                    CHECKIMAGE_TYPE=checkimage_type,
                    CHECKIMAGE_NAME=checkimage_name
                    ))
        
        kwargs = dict(code='SExtractor')
        kwargs['config'] = configdefaults.copy()
        kwargs['temp_path'] = '.'
        kwargs['params'] = self.internals['params']
        
        if config is not None:
            assert isinstance(config,dict)
            kwargs['config'].update(config)
             
        return kwargs
    
    def run_sex(self,catroot,config=None,checks=None,cleanafter=False):
        """ """
        kwargs = self.get_sex_kwargs(catroot,config,checks)
        sextractor = aw.api.Astromatic(**kwargs)
        sextractor.run_frames(self.files['kmage'],frames=[1])
        catname = kwargs['config']['CATALOG_NAME']
        
        if cleanafter:
            self.cleanaftersex(kwargs['config'])
        
        return catname
    
    def cleanaftersex(self,config):
        if 'CHECKIMAGE_TYPE' in config:
            checkimage_name = config['CHECKIMAGE_NAME']
            checkimages = st.split(checkimage_name,',')
            for chimg in checkimages:
                if os.path.exists(chimg):
                    os.system('rm %s' % chimg)
    
    def load_catalog(self,catpath):
        """ """
        catalog = aw.utils.ldac.get_table_from_ldac(catpath)
        return catalog
    
        