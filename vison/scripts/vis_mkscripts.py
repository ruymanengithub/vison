#! $HOME/SOFTWARE/anaconda2/envs/vison/bin/ python
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
import copy
from optparse import OptionParser
import datetime
import string as st
import sys

import vison
from vison.datamodel import scriptic as sc
from vison.campaign import campaign, minicampaign
from vison.support import context
from vison.campaign.chronos import get_test_duration
from vison.support import vjson
# END IMPORT


def f_write_script(struct, filename, outpath, elvis):

    script = sc.Script(structure=struct, elvis=elvis)
    script.build_cargo()
    script.write(filename)

    stdout = os.popen('md5sum %s' % filename).read()
    xsum = st.split(stdout)[0]

    os.system('mv %s %s/' % (filename, outpath))

    return xsum


def scwriter(toWrite, test_generator, outpath, equipment, elvis=context.elvis):
    """ """

    datetag = (datetime.datetime.now()).strftime('%d%b%y')
    versiontag = vison.__version__

    checksumf = 'CHECK_SUMS_%s.txt' % datetag
    inventoryf = 'TESTS_INVENTORY_%s.txt' % datetag
    checksums = []

    if not os.path.exists(outpath):
        os.system('mkdir %s' % outpath)

    test_sequence = test_generator(equipment, toWrite, elvis=elvis)

    Nframes = 0
    duration = 0

    print '\nWRITING SCRIPTS...'
    
    summaries = []
    
    for test in test_sequence.keys():
        print '%s...' % test
        testobj = copy.deepcopy(test_sequence[test])
        structtest = testobj.build_scriptdict(elvis=elvis)
        testduration = get_test_duration(structtest)

        duration += testduration

        iNcols = structtest['Ncols']
        frameslist = [structtest['col%i' % i]['frames']
                      for i in range(1, iNcols+1)]
        iNframes = np.sum(frameslist)

        summaries.append('%s: %i [%.2f min] cols: %s' % (test, iNcols, testduration,
                                                  frameslist.__repr__()))
        Nframes += iNframes

        ffile = 'vis_CalCamp_%s_%s_v%s.xlsx' % (test, datetag, elvis)
        xsum = f_write_script(structtest, ffile, outpath, elvis)
        checksums.append((ffile, xsum))
    
    with open(os.path.join(outpath, inventoryf),'w') as f1:
    
        print >> f1, 'Scripts written on %s' % datetag
        print >> f1, 'checksumf: %s' % checksumf
        print >> f1, 'vison version: %s\n' % versiontag
        
        for item in summaries:
            print >> f1, item
        
        print >> f1, '\n %i Frames Total' % Nframes
        print >> f1, '\n %.2f Minutes Total' % duration


    # WRITING CHECKSUMS

    with open(os.path.join(outpath, checksumf), 'w') as f2:
        for item in checksums:
            print >> f2, '%-60s\t%s' % item


if __name__ == '__main__':
    
    #elvis = '6.5.X'
    #camptype = 'Full'  # Mini/Full
    #outpath = 'MiniCal_data02/Full_Camp_scripts_26JAN18_E6.5.0'
    #equipment = dict(operator='cpf',
    #                 sn_ccd1='D01',
    #                 sn_ccd2='D02',
    #                 sn_ccd3='D03',
    #                 sn_roe='SIM',
    #                 sn_rpsu='PSUs')
    #toWrite_def = OrderedDict(BIAS01=0, DARK01=0, CHINJ00=0, CHINJ01=0, CHINJ02=0,
    #                          FLAT01=0, FLAT02=0, PTC01=0, PTC02WAVE=0, PTC02TEMP=0, NL01=0,
    #                          PSF01=0, PSF02=0,
    #                          TP00=0, TP01=0, TP02=0,
    #                          PERSIST01=0, FOCUS00=0)
    
    #inputs = OrderedDict(camptype=camptype,elvis=elvis,outpath=outpath,
    #                     equipment=equipment,toWrite=toWrite)
    
    
    parser = OptionParser()
    parser.add_option("-j","--json",dest="json",default='',help="json file with inputs")
    
    (options, args) = parser.parse_args()
    
    if options.json == '':
        parser.print_help()
        sys.exit()
    
    inputs = vjson.load_jsonfile(options.json,useyaml=True)
    
    # MISSING: INPUTS VALIDATION
     
    if inputs['camptype'] == 'Full':
        test_generator = campaign.generate_test_sequence
    elif inputs['camptype'] == 'Mini':
        test_generator = minicampaign.generate_reduced_test_sequence
    
    
    scwriter(inputs['toWrite'], test_generator, inputs['outpath'], 
             inputs['equipment'],inputs['elvis'])
