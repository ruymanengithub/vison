"""


Functions and variables related to flags for vison.

Created on Wed Sep 20 17:05:00 2017

:author: Ruyman Azzollini

"""
# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
from collections import OrderedDict
# END IMPORT

flagsdict = OrderedDict(MISSDATA=[2**0],
                        POORQUALDATA=[2**1],
                        BADDATA=[2**2],
                        WRONGTEST=[2**3],
                        WRONG_ROE_CONF=[2**4],
                        SOURCEON=[2**5],
                        SOURCEOFF=[2**6],
                        WRONGSOURCE=[2**7],
                        HK_OOL=[2**8],
                        CCDTEMP_OOL=[2**9],
                        UNDERPERFORM=[2**10],
                        REQNOTMET=[2**11],
                        BADDATAPROD=[2**12],
                        DEADCCD=[2**13],
                        DEADQUAD=[2**14],
                        LOWMEMORY=[2**15],
                        TOOSLOW=[2**16],
                        UNKNOWN=[2**17],
                        RON_OOL=[2**18],
                        FLUENCE_OOL=[2**19],
                        FLUENCEGRAD_OOL=[2**20],
                        FOCUS_OOL=[2**21],
                        BGD_OOL=[2**22],
                        CCDSATUR=[2**23],
                        FLUX_OOL=[2**24],
                        STARS_MISSING=[2**25],
                        SUBTASKCRASH=[2**26],
                        BADFIT=[2**27],
                        SATURATION=[2**28],
                        LOSTPIXELS=[2**29],
                        NONSATURATION=[2**30])


class Flags(object):
    """ """

    def __init__(self, indict=None):
        if indict is None:
            indict = flagsdict
        self.value = 0
        self.dict = indict.copy()

    def isflagon(self, flag):
        return bool(self.value & self.dict[flag][0] == self.dict[flag][0])

    def add(self, flag):
        self.value |= self.dict[flag][0]

    def getFlagsOnList(self):
        return [key for key in self.dict.keys() if self.isflagon(key)]


if __name__ == '__main__':

    flags = Flags()

    keys = flags.dict.keys()
    flagvals = [flags.dict[key][0] for key in keys]
    maxflag = np.max(flagvals)

    maxpow = int(np.log10(maxflag * 1.) / np.log10(2.))

    for key in keys:
        flags.add(key)

    test = 0
    for i in range(maxpow + 1):
        test = test | 2**(i * 1)

    assert test == flags.value

    areOn = True

    for key in keys:
        areOn &= flags.isflagon(key)

    assert areOn

    print 'flags.py: Tests Ran Succesfully'
