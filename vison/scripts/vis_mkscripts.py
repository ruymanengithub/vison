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

import vison
from vison.datamodel import scriptic as sc
#from vison.point import FOCUS00,PSF0X
#from vison.dark import BIAS01,DARK01
#from vison.flat import NL01, PTC0X, FLAT0X
#from vison.inject import CHINJ01,CHINJ02
#from vison.pump import TP01, TP02
#from vison.other import PERSIST01 as PER01
#from vison.point import lib as polib
from vison.pipe import campaign, minicampaign
#from vison.pipe import lib as pilib
from vison.support import context

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


def scwriter(toWrite,test_generator,outpath,equipment,elvis=context.elvis):
    """ """
    
    datetag = (datetime.datetime.now()).strftime('%d%b%y')
    versiontag = vison.__version__
    
    checksumf = 'CHECK_SUMS_%s.txt' % datetag
    inventoryf = 'TESTS_INVENTORY_%s.txt' % datetag
    checksums = []
    
    if not os.path.exists(outpath):
        os.system('mkdir %s' % outpath)
        
    test_sequence = test_generator(equipment,toWrite,elvis=elvis)
    
    Nframes = 0
    
    f1 = open(os.path.join(outpath,inventoryf),'w')
    print >> f1, 'Scripst written on %s' % datetag
    print >> f1, 'checksumf: %s' % checksumf
    print >> f1, 'vison version: %s\n' % versiontag
    
    for test in test_sequence.keys():
        structtest = test_sequence[test]
        
        iNcols = structtest['Ncols']
        frameslist = [structtest['col%i' % i]['frames'] for i in range(1,iNcols+1)]
        iNframes = np.sum(frameslist)
        
        summary = '%s: %i cols: %s' % (test,iNcols,frameslist.__repr__())
        print >> f1, summary
        
        Nframes += iNframes
        
        ffile = 'vis_CalCamp_%s_%s_v%s.xlsx' % (test,datetag,elvis)
        xsum = f_write_script(structtest,ffile,outpath,elvis)
        checksums.append((ffile,xsum))
    
    print >> f1, '\n %i Frames Total' % Nframes
    
    f1.close()
    
    # WRITING CHECKSUMS
            
    f2 = open(os.path.join(outpath,checksumf),'w')
    for item in checksums:
        print >> f2, '%-60s\t%s' % item
    f2.close()



if __name__ =='__main__':
    
    
    elvis = '6.5.X'

    camptype = 'Mini' # Mini/Full    
    
    outpath = 'MiniCal_scripts_19JAN18_E6.5.0'
    
    equipment = dict(operator = 'cpf',
    sn_ccd1 = 'D01',
    sn_ccd2 = 'D02',
    sn_ccd3 = 'D03',
    sn_roe= 'SIM?',
    sn_rpsu = 'PSUs')
    
    toWrite_def = OrderedDict(BIAS01=0,DARK01=0,CHINJ00=0,CHINJ01=0,CHINJ02=0,
                      FLAT01=0,FLAT02=0,PTC01=0,PTC02WAVE=0,PTC02TEMP=0,NL01=0,
                      PSF01=0,PSF02=0,
                      TP00=0,TP01=0,TP02=0,
                      PERSIST01=0,FOCUS00=0)
    
    if camptype == 'Full':
        
        #toWrite = toWrite_def.copy()
        #toWrite.update(OrderedDict(BIAS01=1,DARK01=1,CHINJ01=1,CHINJ02=1,
        #              FLAT01=1,FLAT02=1,PTC01=1,PTC02WAVE=1,PTC02TEMP=1,NL01=1,
        #              PSF01=1,PSF02=1,
        #              TP01=1,TP02=1,
        #              PERSIST01=1,FOCUS00=1))
        
        toWrite = OrderedDict(BIAS01=1,DARK01=1,CHINJ01=1,CHINJ02=1,
                     FLAT01=1,FLAT02=1,PTC01=1,PTC02WAVE=1,PTC02TEMP=1,NL01=1,
                     PSF01=1,PSF02=1,
                     TP01=1,TP02=1,
                     PERSIST01=1,FOCUS00=1)
        
        test_generator = campaign.generate_test_sequence
        
    elif camptype == 'Mini':
        
        toWrite = toWrite_def.copy()
        toWrite.update(OrderedDict(BIAS01=1,DARK01=1,CHINJ00=1,TP00=1,
                      FLAT01=1,FLAT02=1,FOCUS00=1,PSF01=True))
        
        test_generator = minicampaign.generate_reduced_test_sequence
    
    scwriter(toWrite,test_generator,outpath,equipment,elvis)
    
