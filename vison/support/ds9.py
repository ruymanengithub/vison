#! /usr/bin/env python

"""

NEEDS REVISION
TODO: UPDATE

Module with functions related to the handling of ds9 from python
through XPA.


:History:
Created on Thu Mar 17 13:18:10 2016

@author: Ruyman Azzollini

"""

# IMPORT STUFF
import os
import string
import time
import numpy as np
# END IMPORT


class ds9class(object):
    """A very simple class to handle ds9 through xpa."""
    
    def __init__(self):
        self.ego = 'ds9'
        self.MaxOpenTime = 10 # Maximum time to open ds9 in seconds
        self.OpenCheckInterval = 0.2 # seconds
    
    def launch(self):
        """Launches ds9"""

        
        os.system('ds9 &')
        
        startTime = time.time()
        while True:
            time.sleep(self.OpenCheckInterval)
            if self.isOpen():
                return
            if time.time() - startTime > self.MaxOpenTime:
                raise RuntimeError('Could not open ds9 window; timeout')
            
    def isOpen(self):
            """Return True if this ds9 window is open
            and available for communication, False otherwise.
            """
            
            try: 
                shit = self.xpaget('mode')
                return True
            except RuntimeError:
                return False


    def xpaset(self,cmd):
        """Executes xpaset."""
        # IMPORT STUFF
        import os
        # END IMPORT
        execline = 'xpaset -p %s %s' % (self.ego,cmd)
        os.system(execline)

    
    def xpaget(self,cmd):
        """Executes xpaget and retrieves the stdout. If an error happens, an exception is raised."""

        execline = 'xpaget %s %s' % (self.ego,cmd)
        retrieve = os.popen4(execline)
        stdout = retrieve[1].readlines()
        retrieve[0].close()
        retrieve[1].close()
        diderror = string.rfind(stdout[0],'ERROR')
        if diderror > -1:
            raise RuntimeError(stdout[0])
        else :
            return stdout            

    def zoomhere(self,x,y,zoom):
        """Zooms in on given coordinates of display (ds9)."""
        
        try : 
            self.xpaset('zoom to fit')
            self.xpaset('zoom %i' % zoom)
            if (not np.isnan(x)) and (not np.isnan(y)):
	        self.xpaset('pan to %i %i' % (x,y))
        except RuntimeError : return
