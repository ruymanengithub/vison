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
# END IMPORT

flagsdict = {'MISSDATA': [2**0L], 'BADQUALDATA':[2**1L],'HK_OOL':[2**2L],
            'CCDTEMP_OOL':[2**3L],'WRONGTEST':[2**4L],'SOURCEOFF':[2**5L],
            'WRONGSOURCE':[2**6L],'FLUENCE_OOL':[2**7L],
            'WRONG_ROE_CONF':[2**8L],'FOCUS_OOL':[2**9L],
            'RON_OOL':[2**10L]}



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
    
    
if __name__ == '__main__':
    
    flags = Flags()
    
    keys = flags.dict.keys()
    flagvals = [flags.dict[key][0] for key in keys]
    maxflag = np.max(flagvals)
    maxpow = int(np.log10(maxflag)/np.log10(2))
    
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


