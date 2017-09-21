#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

DataDict Class : holds data and results across sub-tasks of a "task" (Test).
This is the core data-structure used to do analysis and report results.


Created on Thu Sep 21 16:47:09 2017

@author: raf
"""


class DataDict():
    """ """
    
    
    def __init__(self,meta=dict()):
        """ """
        
        self.meta = meta
        self.xobs = None
        self.xccd = None
        self.xquad = None
        self.xsubquad = None
        self.dps = dict()
        
    
#ngdd['xobs'] = convert_to_DFrame(DataDict,keys) # indexes: OBSID
    #ngdd['xccd'] = convert_to_MI_DFrame(DataDict,commkeys,CCD=True,Quad=False,Spot=False) # indexes: OBSID, CCD
    #ngdd['xquad'] = None # indexes: OBSID, CCD, Quad
    #ngdd['xsubquad'] = None # indexes: OBSID, CCD, Quad, Spot
    
def useCases():
    """ 
    
    #TODO:
        
        # create a 
    
    """
    
    
    


if __name__ == '__main__':
    
    useCases()