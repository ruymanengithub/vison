#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 16:06:40 2019

@author: raf
"""

# IMPORT STUFF
from pdb import set_trace as stop
from collections import OrderedDict
import os
import copy
import numpy as np
import os
from astropy import table
import glob

from vison.support.files import cPickleRead, cPickleDumpDictionary
from vison.support import vjson
from vison.datamodel import core as vcore
from vison.support import vcal
from vison.fpa import fpa as fpamod

# END IMPORT

vcalpath = 'FPA_CAL'
vcalfile = os.path.join(vcalpath,'CCD_CAL_CONVERSIONS_ALL_BLOCKS.xlsx')

class MetaCal(object):
    
    def __init__(self, respathroot, outpathroot):
        """ """
        
        self.blocks = fpamod.all_blocks
        self.flight_blocks = fpamod.flight_blocks
        self.CCDs = [1,2,3]
        self.Quads = ['E','F','G','H']
        self.respathroot = respathroot
        self.outpath = outpathroot
        self.inventory = OrderedDict()
        self.results = OrderedDict()
        self.ParsedTable = None
        self.roeVCals = dict()        
        self.init_roeVCals()
        
    def init_roeVCals(self):
        
        for block in self.blocks:
            
            ROE_SN = fpamod.ROE_SNs[block]
            
            roevcal = vcal.RoeVCal(vcalfile)
            try: roevcal.load_VCal(ROE_SN)
            except IOError:
                continue
                
            self.roeVCals[block] = copy.deepcopy(roevcal)
            

    def _update_inventory(self, block, testline):
    
        test = testline[0]
        session = testline[1]
        repeat = testline[2]
        
        sesspath = os.path.join(self.respathroot,block,'ANALYSIS','results_atCALDATA',
                               session)
        
        alltestpaths = glob.glob(os.path.join(sesspath,'%s*' % test))
        
        if len(alltestpaths)>1:
            testpath = '%s.%i' % (test,repeat)
        else:
            testpath = test
            
        respath = os.path.join(sesspath, testpath)
        
        DDfile = os.path.join(respath,'%s_DataDict.pick' % test)
        
        
        try:
            assert os.path.exists(respath)
        except:
            stop()
        try:
            assert os.path.exists(DDfile)
        except:
            stop()
        
        if block not in self.inventory:
            self.inventory[block] = OrderedDict()
        if test not in self.inventory[block]:
            self.inventory[block][test] = []
        
        self.inventory[block][test].append(OrderedDict())
        
        ii = self.inventory[block][test][-1]
        
        ii['DD'] = DDfile
        dd = cPickleRead(DDfile)
        ii['dd'] = copy.deepcopy(dd)
        ii['resroot'] = respath
        ii['session'] = session
        ii['repeat'] = repeat
        
        return None
    
    
    def load_block_results(self,inventoryfile):
        
        inventraw = vjson.load_jsonfile(inventoryfile, useyaml=True)    
        
        
        for block in self.blocks:
            if block in inventraw['inventory'].keys():
                rawcargo = inventraw['inventory'][block]
                
                for testline in rawcargo:                
                    self._update_inventory(block, testline)
    
        return None
    
    def parse_block_results(self,*args,**kwargs):
        raise NotImplementedError("Subclass must implement abstract method")
    
    def dump_aggregated_results(self):
        raise NotImplementedError("Subclass must implement abstract method")
    
    
    def stack_dd(self, dd, cols, index2stack, indices2keep, stacker='median'):
        """ """
        
        #indices = copy.deepcopy(dd.indices)
        
        stkdd = vcore.DataDict()
        
        
        stacker_functions = dict(median=np.median,
                                 mean=np.mean)
        
        fstacker = stacker_functions[stacker]
        
        for ixc, icol in enumerate(cols):
            
            idtype = dd.mx[icol][:].dtype
            
            iindices = dd.mx[icol].indices
            
                            
            niindices = []
                                        
            niindices = [index for index in iindices \
                         if index.name in indices2keep]
            
            # STACKING
            
            _ix2stack = iindices.names.index(index2stack)
            
            if idtype.char not in ['S', 'O']:
                
                imx = fstacker(dd.mx[icol][:],axis=_ix2stack,keepdims=True)
            
            else:
                slicer = []
                for _ix, _i in enumerate(iindices):
                    if _i.name == index2stack:
                        slicer.append([0])
                    else:
                        slicer.append(slice(_i.len))
                imx = dd.mx[icol][tuple(slicer)]
            
           
            # TRIMMING UNNECESSARY AXES
                       
            
            if idtype.char not in ['S', 'O']:
                
                _ix2ax = [_ix for _ix,ix in enumerate(iindices) if \
                          ix.name not in indices2keep]
                
                imx = fstacker(imx,axis=_ix2ax,keepdims=False)
            else:
                
                slicer = []
                
                for _ix, _i in enumerate(iindices):
                    if _i.name not in indices2keep and _i.name != index2stack:
                        slicer.append(0)
                
                imx = imx[tuple(slicer)]
            
            # INDEXING
            
            niindices = []
            
            for ix in iindices:
                if ix.name not in indices2keep:
                    pass
                else:
                    if ix.name == index2stack:
                        _ix = copy.deepcopy(ix)
                        _ix.len=1
                        _ix.vals=[ix.vals[0]]
                        niindices.append(_ix)
                    else:
                        niindices.append(copy.deepcopy(ix))
            
            stkdd.addColumn(imx,name=icol,indices=niindices)

            
        return stkdd
    
    def stackTables(self, t1, t2):
        """ """        
        return table.vstack([t1,t2])

        
    def get_FPAMAP_from_PT(self,PT,extractor):
        """ """
        
        M = OrderedDict()
        
        NSLICES = fpamod.NSLICES
        NCOLS = fpamod.NCOLS
        
        for jY in range(NSLICES):
            for iX in range(NCOLS):
                Ckey  = 'C_%i%i' % (jY+1,iX+1)
                M[Ckey] = OrderedDict()
                
                locator = fpamod.FPA_MAP[Ckey]
                block = locator[0]
                CCDk = locator[1]
                
                for Q in self.Quads:
                    
                    M[Ckey][Q] = extractor(PT,block,CCDk,Q)
        
        return M
        
        
        