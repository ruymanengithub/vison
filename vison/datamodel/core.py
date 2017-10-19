#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

DataDict Class : holds data and results across sub-tasks of a "task" (Test).
This is the core data-structure used to do analysis and report results.


Created on Thu Sep 21 16:47:09 2017

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
from collections import OrderedDict

from vison.pipe import lib as pilib
# END IMPORT



class Index(object):
    """ """
    
    def __init__(self,name,vals=[],N=0):
        """ """
        self.name = name
        if len(vals) !=0:
            self.vals = vals
            self.len = len(self.vals)
        elif (len(vals) == 0) and N!=0:
            self.vals = np.arange(N).tolist()
            self.len = N
        else:
            raise TypeError
    def __str__(self):
        return '"%s": %i' % (self.name,self.len)




class Column(object):
    """ """
    
    def __init__(self,array,name,indices):
        """ """
        
        self.array = array.copy()
        self.shape = self.array.shape
        assert isinstance(indices,list)
        for index in indices:
            assert isinstance(index,Index)
            
        assert len(self.shape) == len(indices)
        for i in range(len(self.shape)):
            assert self.shape[i] == indices[i].len, stop()
                
        self.indices = indices 
        self.name = name
    
    
    def __call__(self):
        return self.array
    
    def __str__(self):
        return '%s: %s' % (self.name,self.array.__str__())
    
    #def __str__(self):
    #    return self.name

class DataDict():
    """ """
    
    
    def __init__(self,meta=dict()):
        """ """
        
        self.meta = meta
        self.mx = OrderedDict()
        self.colnames = []
        self.indices = []
        self.products = dict() # data products
        
    def loadExpLog(self,explog):
        """ """
        
        CCDs = [1,2,3]
        
        ObsID = explog['ObsID'].data
        uObsID = np.unique(ObsID)
        Nobs = len(uObsID)
        
        ccdcol = explog['CCD']
        
        ObsIndex = Index('ix',N=Nobs)
        commIndices = [ObsIndex,Index('CCD',CCDs)]
        
        self.addColumn(uObsID,'ObsID',[ObsIndex])
        
        for key in explog.colnames:
                        
            if key == 'ObsID':                
                continue
                
            arrlist = []
            for CCDindex in CCDs:
                CCDkey = 'CCD%i' % CCDindex
                arrlist.append(explog[key].data[np.where(ccdcol==CCDkey)])
            array = np.array(arrlist).transpose() 
            self.addColumn(array,key,commIndices)
                    
        return None
        

    
    def addColumn(self,array,name,indices):
        """ """
        
        column = Column(array,name,indices)
        
        self.mx[column.name] = column
        self.colnames.append(column.name)
        colindices = column.indices
        
        selfindnames = [index.name for index in self.indices]
        
        for ic,index in enumerate(colindices):
            if index.name in selfindnames:
                ix = selfindnames.index(index.name)
                assert np.all(index.vals == self.indices[ix].vals)
                assert ic == ix
            else:
                self.indices.append(index)
        
        
    def dropColumn(self,):
        """ """
        
        
    
    def saveToFile(self,):
        """ """
        
    
    
    
def useCases():
    """ 
    
    #TODO:
        
        # create a DataDict object from an exposure log.
        # add a column indexed by ObsID, CCD and Quad
        # drop a column
        # save to a text / excel file
        # create a column from an operation on several columns with different dimensions
        
    
    """
    
    dpath = '/home/raf/WORK/EUCLID/CALIBRATION/PipelineDevel/TEST_DATA/24_Feb_80/'
    explogf = os.path.join(dpath,'EXP_LOG_240280.txt')
    
    explog = pilib.loadexplogs(explogf,elvis='6.3.0',addpedigree=False,datapath=None)
    
    dd = DataDict()
    
    dd.loadExpLog(explog)
    
    
    stop()
    


if __name__ == '__main__':
    
    useCases()