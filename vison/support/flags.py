"""

VIS Ground Calibration


Functions and variables related to flags for vison.

Created on Wed Sep 20 17:05:00 2017

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""
# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
from collections import OrderedDict
# END IMPORT

flagsdict = OrderedDict(MISSDATA=[2**0L],
            POORQUALDATA=[2**1L],
            BADDATA=[2**2L],
            WRONGTEST=[2**3L],
            WRONG_ROE_CONF=[2**4L], 
            SOURCEON=[2**5L],            
            SOURCEOFF=[2**6L],
            WRONGSOURCE=[2**7L],
            FLUENCE_OOL=[2**8L],
            FOCUS_OOL=[2**9L],
            BGDHI=[2**10L],
            CCDSATUR=[2**11L],
            HK_OOL=[2**12L],
            CCDTEMP_OOL=[2**13L],
            RON_OOL=[2**14L],
            UNDERPERFORM=[2**15L],
            REQNOTMET=[2**16L],               
            BADDATAPROD=[2**17L],
            DEADCCD=[2**18L],
            DEADQUAD=[2**19L],
            LOWMEMORY=[2**20L],
            TOOSLOW=[2**21L],
            UNKNOWN=[2**22L])



class Flags(object):
    """ """
    
    def __init__(self,):
        self.value = 0L
        self.dict = flagsdict
    
    def isflagon(self,flag):
        if self.value & self.dict[flag][0] == self.dict[flag][0]: 
            return True
        else: return False
        
    def add(self,flag):
        self.value |= self.dict[flag][0]
    
    def getFlagsOnList(self):
        return [key for key in self.dict.keys() if self.isflagon(key)]
        
    
    
if __name__ == '__main__':
    
    flags = Flags()
    
    keys = flags.dict.keys()
    flagvals = [flags.dict[key][0] for key in keys]
    maxflag = np.max(flagvals)
    
    maxpow = int(np.log10(maxflag*1.)/np.log10(2.))
    
    for key in keys:
        flags.add(key)
    
    test=0
    for i in range(maxpow+1):
        test = test | 2**(i*1L)
        
    
    assert test == flags.value
    
    areOn = True
    
    for key in keys:
        areOn &= flags.isflagon(key)
    
    assert areOn == True
    
    
    print 'flags.py: Tests Ran Succesfully'


