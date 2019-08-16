#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 11:28:28 2019

@author: raf
"""

# IMPORT STUFF


# END IMPORT


flight_blocks = ['BORN','CURIE','DIRAC','ERWIN','FOWLER','GUYE','SKLODOWSKA',
          'JULES2','KRAMERS','LORENTZ','MAX','NIELS']
spare_blocks = ['OWEN','EINSTEIN']

all_blocks = flight_blocks + spare_blocks


NSLICES = 6
NCOLS = 6

FPA_MAP = dict(
# SLICE_1
        C_11 = ('BORN','CCD1',(1,0)),
        C_12 = ('BORN','CCD2',(1,0)),
        C_13 = ('BORN','CCD3',(1,0)),
        C_14 = ('DIRAC','CCD3',(0,1)),
        C_15 = ('DIRAC','CCD2',(0,1)),
        C_16 = ('DIRAC','CCD1',(0,1)),
# SLICE_2
        C_21 = ('GUYE','CCD1',(1,0)),
        C_22 = ('GUYE','CCD2',(1,0)),
        C_23 = ('GUYE','CCD3',(1,0)),
        C_24 = ('FOWLER','CCD3',(0,1)),
        C_25 = ('FOWLER','CCD2',(0,1)),
        C_26 = ('FOWLER','CCD1',(0,1)),
# SLICE_3 
        C_31 = ('ERWIN','CCD1',(1,0)),
        C_32 = ('ERWIN','CCD2',(1,0)),
        C_33 = ('ERWIN','CCD3',(1,0)),
        C_34 = ('KRAMERS','CCD3',(0,1)),
        C_35 = ('KRAMERS','CCD2',(0,1)),
        C_36 = ('KRAMERS','CCD1',(0,1)),
# SLICE_4 
        C_41 = ('JULES2','CCD1',(1,0)),
        C_42 = ('JULES2','CCD2',(1,0)),
        C_43 = ('JULES2','CCD3',(1,0)),
        C_44 = ('SKLODOWSKA','CCD3',(0,1)),
        C_45 = ('SKLODOWSKA','CCD2',(0,1)),
        C_46 = ('SKLODOWSKA','CCD1',(0,1)),
# SLICE_5 
        C_51 = ('LORENTZ','CCD1',(1,0)),
        C_52 = ('LORENTZ','CCD2',(1,0)),
        C_53 = ('LORENTZ','CCD3',(1,0)),
        C_54 = ('MAX','CCD3',(0,1)),
        C_55 = ('MAX','CCD2',(0,1)),
        C_56 = ('MAX','CCD1',(0,1)),
# SLICE_6
        C_61 = ('NIELS','CCD1',(1,1)),
        C_62 = ('NIELS','CCD2',(1,1)),
        C_63 = ('NIELS','CCD3',(1,1)),
        C_64 = ('CURIE','CCD3',(0,1)),
        C_65 = ('CURIE','CCD2',(0,1)),
        C_66 = ('CURIE','CCD1',(0,1)),
    )

BLOCK_SNs = dict(BORN=1,
            CURIE=4,
            DIRAC=3,
            ERWIN=2,
            FOWLER=5,
            GUYE=6,
            HEISENBERG=7,
            SKLODOWSKA=7,
            JULES=8,            
            JULES2=8,
            KRAMERS=9,
            LORENTZ=12,
            MAX=10,
            NIELS=11,
            OWEN=13,
            EINSTEIN=14)

ROE_SNs = dict(BORN=1,
            CURIE=4,
            DIRAC=5,
            ERWIN=2,
            FOWLER=6,
            GUYE=3,
            HEISENBERG=7,
            SKLODOWSKA=7,
            JULES=8,            
            JULES2=8,
            KRAMERS=9,
            LORENTZ=12,
            MAX=10,
            NIELS=11,
            OWEN=13,
            EINSTEIN=14)

RPSU_SNs = dict(BORN=1,
            CURIE=4,
            DIRAC=5,
            ERWIN=2,
            FOWLER=6,
            GUYE=3,
            HEISENBERG=7,
            SKLODOWSKA=7,
            JULES=8,            
            JULES2=8,
            KRAMERS=9,
            LORENTZ=12,
            MAX=10,
            NIELS=11,
            OWEN=13,
            EINSTEIN=14)


def flip_img(img,flip):
    
    if flip[0] == 1:
        img = img[::-1,:].copy()
    if flip[1] == 1:
        img = img[:,::-1].copy()
    return img


class FPA(object):
    """ """
        
    def __init__(self):
        """ """
        self.NSLICES=6
        self.NCOLS=6
        self.flight_blocks = flight_blocks
        self.spare_blocks = spare_blocks
        self.all_blocks = all_blocks        
        self.placeholderdata = None
        self.FPA_MAP = FPA_MAP
    
    def flip_img(self,img,flip):
        return flip_img(img,flip)
    
    def get_Ckey_from_BlockCCD(self, block, CCD):
        """ """
        CCDk = 'CCD%i' % CCD
        
        for jY in range(1,self.NSLICES+1):
            for iX in range(1,self.NCOLS+1):
                Ckey = 'C_%i%i' % (jY,iX)
                locator = self.FPA_MAP[Ckey]
                if (block == locator[0]) and (CCDk == locator[1]):
                    return Ckey
        return None
    
    
    def assign_ccd_data(self,blocksdata,Ckey):
        """ """        
        block, CCDk, _ = self.FPA_MAP[Ckey]
        try:
            return blocksdata[block][CCDk].copy()
        except KeyError:
            return self.placeholderdata
            
    
    def assign_ccd_img(self,blocksdata,Ckey):
        """ """        
        block, CCDk, flip = self.FPA_MAP[Ckey]
        try:
            return self.flip_img(blocksdata[block][CCDk],flip)
        except KeyError:
            return self.placeholderdata
    
    
    def build_fpa_data(self,blocksdata,kind):
        """ """
        
        fpadata = dict()
        
        if kind == 'data':
            _method = self.assign_ccd_data
        elif kind == 'image':
            _method = self.assign_ccd_img
    
        for jY in self.NSLICES:
            for iX in self.NCOLS:                
                Ckey = 'C_%i%i' % (jY,iX)                
                fpadata[Ckey] = _method(blocksdata,Ckey)
        
        return fpadata
        
        