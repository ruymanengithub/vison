"""Functions and variables related to flags for moments.py"""

import numpy as num

allflags = {'NOPETRO': 2**0L, 'MANYPETRO':2**1L, 'NEGPETRO':2**2L,\
'BADRADIAL':2**3L,'MAXITER_AS':2**4L,'NOCONC':2**5L,'MANYCONC':2**6L,\
'NORADIAL':2**8L,'USEPETROMSK':2**9L,'PETROMSKTRUNC':2**10L,\
'MAXITER_AXAS':2**11L,'NOM20':2**12L,'NEGGROWTH':2**13L,'BLANK':2**14L,\
'NONCHECKEDRADIAL':2**15L,'NOFFACTOR':2**16L,'NOCLUMPS':2**17L,\
'NOPEAKS':2**18L,'SHORTSKYRAD':2**19L}

#for key in allflags.keys(): allflags[key] = num.long(allflags[key])

def isflagon(allflags,flag):
    if allflags & flag == flag : return True
    else : return False

def addflag(allflags,flag):
    """Adds a flag to a flag variable."""
    # IMPORT STUFF
    from pdb import set_trace as stop
    # END IMPORT
    #if type(allflags) != type(flag) : stop()
    allflags = allflags | flag    
    return allflags
    
