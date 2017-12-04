"""

VIS Ground Calibration Campaign


Automatically Generating Calibration Campaign Scripts.

Created on Fri Sep 08 12:03:00 2017

:autor: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import os
from collections import OrderedDict

from vison.datamodel import elvis as elv
from vison.datamodel import scriptic as sc
#from vison.point import FOCUS00,PSF0X
#from vison.dark import BIAS01,DARK01
#from vison.flat import NL01, PTC0X, FLAT0X
#from vison.inject import CHINJ01,CHINJ02
#from vison.pump import TP01, TP02
#from vison.other import PERSIST01 as PER01
#from vison.point import lib as polib
from vison.pipe import campaign

from vison.ogse.ogse import tFWC_flat,tFWC_point

import datetime
import string as st

# END IMPORT

def f_write_script(struct,filename,outpath,elvis): 
    
    script = sc.Script(structure=struct,elvis=elvis)
    script.build_cargo()
    script.write(filename)
    
    stdout = os.popen('md5sum %s' % filename).read()
    xsum = st.split(stdout)[0]
    
    os.system('mv %s %s/' % (filename,outpath))   
    
    return xsum


def scwriter(toWrite,outpath,equipment,elvis='6.3.0'):
    """ """
    
    
    datetag = (datetime.datetime.now()).strftime('%d%b%y')
    
    checksumf = 'CHECK_SUMS.txt'
    checksums = []
    
    if not os.path.exists(outpath):
        os.system('mkdir %s' % outpath)
        
    test_sequence = campaign.generate_test_sequence(equipment,toWrite,elvis=elvis)    

    
    for test in test_sequence.keys():
        structtest = test_sequence[test]
        ffile = 'vis_CalCamp_%s_%s_v%s.xlsx' % (test,datetag,elvis)
        xsum = f_write_script(structtest,ffile,outpath,elvis)
        checksums.append((ffile,xsum))
    
    
    # WRITING CHECKSUMS
            
    f = open(os.path.join(outpath,checksumf),'w')
    for item in checksums:
        print >> f, '%-60s\t%s' % item
    f.close()



if __name__ =='__main__':
    
    
    elvis = '6.3.0'
    
    
    outpath = 'CAL_scripts_02DEC'
    
    equipment = dict(operator = 'raf',
    sn_ccd1 = 'CCD1TEST',
    sn_ccd2 = 'CCD2TEST',
    sn_ccd3 = 'CCD3TEST',
    sn_roe= 'ROETEST',
    sn_rpsu = 'RPSUTEST')
    
    toWrite = OrderedDict(BIAS01=1,DARK01=1,CHINJ01=1,CHINJ02=1,
                      FLAT01=1,FLAT02=1,PTC01=1,PTC02WAVE=1,PTC02TEMP=1,NL01=1,
                      PSF01=1,PSF02=1,
                      TP01=1,TP02=1,
                      PERSIST01=1,FOCUS00=1)
    
    scwriter(toWrite,outpath,equipment,elvis)
    
