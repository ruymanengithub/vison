#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 11:28:28 2019

@author: raf
"""

# IMPORT STUFF
from pdb import set_trace as stop

# END IMPORT


flight_blocks = ['BORN','CURIE','DIRAC','ERWIN','FOWLER','GUYE','SKLODOWSKA',
          'JULES2','KRAMERS','LORENTZ','MAX','NIELS']
spare_blocks = ['OWEN','EINSTEIN']

all_blocks = flight_blocks + spare_blocks

def get_block_lists(design='initial'):
    
    block_lists = dict()
    
    if design == 'initial':
        
        block_lists['flight_blocks'] = ['BORN','CURIE','DIRAC','ERWIN','FOWLER','GUYE','SKLODOWSKA',
          'JULES2','KRAMERS','LORENTZ','MAX','NIELS']
        block_lists['spare_blocks'] = ['OWEN', 'EINSTEIN']
        block_lists['discarded_blocks'] = ['HEISENBERG','JULES']
        

    elif design == 'final':
        
        block_lists['flight_blocks'] = ['BORN','CURIE','DIRAC','ERWIN','FOWLER','GUYE',
          'JULES2','KRAMERS','LORENTZ','NIELS', 'OWEN', 'EINSTEIN']
        block_lists['spare_blocks'] = ['SKLODOWSKA','MAX']
        block_lists['discarded_blocks'] = ['HEISENBERG','JULES']
        
    block_lists['all_blocks'] = block_lists['flight_blocks'] + block_lists['spare_blocks'] +\
               block_lists['discarded_blocks']
    
    
    
    return block_lists

NSLICES = 6
NCOLS = 6

FPA_MAP_initial = dict(
# SLICE_1
        C_11 = ('BORN','CCD1',(1,0),'15081-08-02'),
        C_12 = ('BORN','CCD2',(1,0),'15081-16-02'),
        C_13 = ('BORN','CCD3',(1,0),'15081-22-02'),
        C_14 = ('DIRAC','CCD3',(0,1),'14313-02-01'),
        C_15 = ('DIRAC','CCD2',(0,1),'14471-01-02'),
        C_16 = ('DIRAC','CCD1',(0,1),'15081-07-02'),
# SLICE_2
        C_21 = ('GUYE','CCD1',(1,0),'14462-13-02'),
        C_22 = ('GUYE','CCD2',(1,0),'14313-17-01'),
        C_23 = ('GUYE','CCD3',(1,0),'14313-23-01'),
        C_24 = ('FOWLER','CCD3',(0,1),'14313-13-01'),
        C_25 = ('FOWLER','CCD2',(0,1),'14311-05-02'),
        C_26 = ('FOWLER','CCD1',(0,1),'14311-20-02'),
# SLICE_3 
        C_31 = ('ERWIN','CCD1',(1,0),'15053-01-02'),
        C_32 = ('ERWIN','CCD2',(1,0),'14352-04-02'),
        C_33 = ('ERWIN','CCD3',(1,0),'14471-16-01'),
        C_34 = ('KRAMERS','CCD3',(0,1),'14471-15-01'),
        C_35 = ('KRAMERS','CCD2',(0,1),'14462-16-01'),
        C_36 = ('KRAMERS','CCD1',(0,1),'14362-20-02'),
# SLICE_4 
        C_41 = ('JULES2','CCD1',(1,0),'14462-02-01'),
        C_42 = ('JULES2','CCD2',(1,0),'14462-10-01'),
        C_43 = ('JULES2','CCD3',(1,0),'14313-06-01'),
        C_44 = ('SKLODOWSKA','CCD3',(0,1),'14313-24-02'),
        C_45 = ('SKLODOWSKA','CCD2',(0,1),'14313-08-01'),
        C_46 = ('SKLODOWSKA','CCD1',(0,1),'14471-20-02'),
# SLICE_5 
        C_51 = ('LORENTZ','CCD1',(1,0),'15053-02-02'),
        C_52 = ('LORENTZ','CCD2',(1,0),'15053-20-02'),
        C_53 = ('LORENTZ','CCD3',(1,0),'15081-22-01'),
        C_54 = ('MAX','CCD3',(0,1),'15081-15-01'),
        C_55 = ('MAX','CCD2',(0,1),'15081-17-01'),
        C_56 = ('MAX','CCD1',(0,1),'15081-21-01'),
# SLICE_6
        C_61 = ('NIELS','CCD1',(1,1),'14471-01-01'),
        C_62 = ('NIELS','CCD2',(1,1),'15081-20-02'),
        C_63 = ('NIELS','CCD3',(1,1),'15053-01-01'),
        C_64 = ('CURIE','CCD3',(0,1),'14462-08-01'),
        C_65 = ('CURIE','CCD2',(0,1),'15081-21-02'),
        C_66 = ('CURIE','CCD1',(0,1),'14301-24-01'),
    )


FPA_MAP_final = dict(
# SLICE_1
        C_11 = ('BORN','CCD1',(1,0),'15081-08-02'),
        C_12 = ('BORN','CCD2',(1,0),'15081-16-02'),
        C_13 = ('BORN','CCD3',(1,0),'15081-22-02'),
        C_14 = ('DIRAC','CCD3',(0,1),'14313-02-01'),
        C_15 = ('DIRAC','CCD2',(0,1),'14471-01-02'),
        C_16 = ('DIRAC','CCD1',(0,1),'15081-07-02'),
# SLICE_2
        C_21 = ('GUYE','CCD1',(1,0),'14462-13-02'),
        C_22 = ('GUYE','CCD2',(1,0),'14313-17-01'),
        C_23 = ('GUYE','CCD3',(1,0),'14313-23-01'),
        C_24 = ('FOWLER','CCD3',(0,1),'14313-13-01'),
        C_25 = ('FOWLER','CCD2',(0,1),'14311-05-02'),
        C_26 = ('FOWLER','CCD1',(0,1),'14311-20-02'),
# SLICE_3 
        C_31 = ('ERWIN','CCD1',(1,0),'15053-01-02'),
        C_32 = ('ERWIN','CCD2',(1,0),'14352-04-02'),
        C_33 = ('ERWIN','CCD3',(1,0),'14471-16-01'),
        C_34 = ('KRAMERS','CCD3',(0,1),'14471-15-01'),
        C_35 = ('KRAMERS','CCD2',(0,1),'14462-16-01'),
        C_36 = ('KRAMERS','CCD1',(0,1),'14362-20-02'),
# SLICE_4 
        C_41 = ('JULES2','CCD1',(1,0),'14462-02-01'),
        C_42 = ('JULES2','CCD2',(1,0),'14462-10-01'),
        C_43 = ('JULES2','CCD3',(1,0),'14313-06-01'),
        C_44 = ('OWEN','CCD3',(0,1),'14362-12-01'),
        C_45 = ('OWEN','CCD2',(0,1),'14311-16-02'),
        C_46 = ('OWEN','CCD1',(0,1),'14462-02-02'),
# SLICE_5 
        C_51 = ('LORENTZ','CCD1',(1,0),'15053-02-02'),
        C_52 = ('LORENTZ','CCD2',(1,0),'15053-20-02'),
        C_53 = ('LORENTZ','CCD3',(1,0),'15081-22-01'),
        C_54 = ('EINSTEIN','CCD3',(0,1),'14471-10-02'),
        C_55 = ('EINSTEIN','CCD2',(0,1),'15081-15-02'),
        C_56 = ('EINSTEIN','CCD1',(0,1),'14471-19-01'),
# SLICE_6
        C_61 = ('NIELS','CCD1',(1,1),'14471-01-01'),
        C_62 = ('NIELS','CCD2',(1,1),'15081-20-02'),
        C_63 = ('NIELS','CCD3',(1,1),'15053-01-01'),
        C_64 = ('CURIE','CCD3',(0,1),'14462-08-01'),
        C_65 = ('CURIE','CCD2',(0,1),'15081-21-02'),
        C_66 = ('CURIE','CCD1',(0,1),'14301-24-01'),
    )

def get_FPA_MAP(design='initial'):
    
    if design == 'initial':
        return FPA_MAP_initial.copy()
    elif design == 'final':
        return FPA_MAP_final.copy()

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
        
    def __init__(self,design='initial'):
        """ """
        self.NSLICES=6
        self.NCOLS=6
        
        block_lists = get_block_lists(design)
        
        self.flight_blocks = block_lists['flight_blocks']
        self.spare_blocks = block_lists['spare_blocks']
        self.all_blocks = block_lists['all_blocks']
        self.placeholderdata = None
        self.FPA_MAP = get_FPA_MAP(design)
    
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
        
        