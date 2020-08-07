#! $HOME/SOFTWARE/anaconda2/envs/vison/bin/ python
"""

Automatically Generating Calibration Campaign Data Acquisition Scripts.
Aimed at ELVIS.

:History:
Created on Fri Sep 08 12:03:00 2017

:autor: Ruyman Azzollini

"""

# IMPORT STUFF
from __future__ import print_function
from pdb import set_trace as stop
import numpy as np
import os
from collections import OrderedDict
import copy
from optparse import OptionParser
import datetime
import string as st
import sys
import pandas as pd

import vison
from vison.datamodel import scriptic as sc
from vison.campaign import campaign, minicampaign, mot, devel
from vison.support import context
from vison.campaign.chronos import get_test_duration
from vison.support import vjson
from vison.support import excel
# END IMPORT


def f_write_script(struct, filename, outpath, elvis):
    """ 

    Function that writes the test structure (as a dictionary) to an
    .xlsx file.

    Calls vison.datamodel.scriptic.Script class.

    """

    script = sc.Script(structure=struct, elvis=elvis)
    script.build_cargo()
    script.write(filename)

    stdout = os.popen('md5sum %s' % filename).read()
    xsum = st.split(stdout)[0]

    os.system('mv %s %s/' % (filename, outpath))

    return xsum


def write_summary_astextfile(summarylines, inventoryfile, meta):
    """Writes a summary of the scripts generated, in .txt format."""

    fulldatetag = meta['date']
    checksumf = meta['checksumfile']
    versiontag = meta['version']
    #outpath = meta['path']
    Nframes = meta['Nframes']
    duration = meta['Total_duration']

    with open(inventoryfile, 'w') as f1:

        print('Scripts written on %s' % fulldatetag, file=f1)
        print('checksumf: %s' % checksumf, file=f1)
        print('vison version: %s\n' % versiontag, file=f1)

        for item in summarylines:
            print(item, file=f1)

        print('\n %i Frames Total' % Nframes, file=f1)
        print('\n %.2f Minutes Total' % duration, file=f1)


def write_summary_asexcel(summarydict, excelf, meta):
    """Writes a summary of the tests generated in .xlsx format.
    Uses vison.support.excel.ReportXL class and Pandas.

    """

    report = excel.ReportXL()
    report.wb.create_sheet('Meta', 0)
    report.dict_to_sheet(meta, 'Meta')
    report.wb.create_sheet('Tests', 1)

    summarydf = pd.DataFrame.from_dict(summarydict)

    report.df_to_sheet(summarydf, 'Tests')
    report.adjust_columns_width('Tests', minwidth=15)
    report.save(excelf)


def scwriter(toWrite, test_generator, outpath, equipment, elvis=context.elvis, CHAMBER=None):
    """Parent function writing the test scripts.

        * calls test_generator to create a sequence of test objects, instantiated with 
          their inputs.
        * the method build_scriptdict() of each test object is called to generate the 
          dictionary-based structure of the test
        * test acquisition durations are estimated
        * the function f_write_script is called to write the script to an excel.
        * a txt file is written with a summary of the tests written.
        * a .xlsx file is written with the same summary.
        * the checksums of the scripts are written to a .txt file.

    :param toWrite: a dictionary with pairs of test_name: flat, where flag is either 1 (create 
        script), or 0 (do not create).
    :type toWrite: dict()
    :param test_generator: a function that instanstiates test objects for a given campaign type.
    :type test_generator: function
    :param outpath: outputs directory. Created on the fly if not already-existent.
    :type outpath: str
    :param equipment: dictionary with equipment serials.
    :type equipment: dict()
    :param elvis: elvis version to write scripts for.
    :type elvis: str
    :param CHAMBER: chamber where the tests will be executed (affects exposure times, for exapmle).
    :type CHAMBER: str

    :return: None

    """

    datetag = (datetime.datetime.now()).strftime('%d%b%y')
    fulldatetag = (datetime.datetime.now()).strftime('%d%b%y %H:%M:%S')
    versiontag = vison.__version__

    if not os.path.exists(outpath):
        os.system('mkdir %s' % outpath)

    checksumf = 'CHECK_SUMS_%s.txt' % datetag
    inventoryf = os.path.join(outpath, 'TESTS_INVENTORY_%s.txt' % datetag)
    excelf = os.path.join(outpath, 'TESTS_SCHEDULER_%s.xlsx' % datetag)

    checksums = []

    test_sequence = test_generator(
        equipment, toWrite, elvis=elvis, CHAMBER=CHAMBER,
        purpose='scripts')

    Nframes = 0
    duration = 0

    print('\nWRITING SCRIPTS...')

    summarylines = []
    summarydict = OrderedDict()
    summarydict['Test'] = []
    summarydict['Ncols'] = []
    summarydict['Nframes'] = []
    summarydict['Duration_min'] = []

    for test in test_sequence.keys():

        print('%s...' % test)
        testobj = copy.deepcopy(test_sequence[test])
        try:
            structtest = testobj.build_scriptdict(elvis=elvis)
        except BaseException:
            stop()

        testduration = get_test_duration(structtest)

        duration += testduration

        iNcols = structtest['Ncols']

        frameslist = [structtest['col%03i' % i]['frames']
                      for i in range(1, iNcols + 1)]
        iNframes = np.sum(frameslist)

        summarylines.append('%s: %i [%.2f min] cols: %s' % (test, iNcols, testduration,
                                                            frameslist.__repr__()))
        Nframes += iNframes

        ffile = 'vis_CalCamp_%s_%s_v%s.xlsx' % (test, datetag, elvis)
        xsum = f_write_script(structtest, ffile, outpath, elvis)
        checksums.append((ffile, xsum))

        summarydict['Test'].append(test)
        summarydict['Ncols'].append(iNcols)
        summarydict['Nframes'].append(iNframes)
        summarydict['Duration_min'].append('%.2f' % testduration)

    metatxt = dict(date=fulldatetag,
                   checksumfile=checksumf,
                   version=versiontag,
                   path=outpath,
                   Nframes=Nframes,
                   Total_duration=duration)

    # WRITING TEXT SUMMARY FILE

    write_summary_astextfile(summarylines, inventoryf, metatxt)

    # WRITING EXCEL SUMMARY FILE

    metaexcel_inputs = [('date', fulldatetag),
                        ('checksumfile', checksumf),
                        ('version', versiontag),
                        ('path', outpath),
                        ('Nframes', Nframes),
                        ('Total_duration', duration)]
    metaexcel = OrderedDict(metaexcel_inputs)

    write_summary_asexcel(summarydict, excelf, metaexcel)

    # WRITING CHECKSUMS

    with open(os.path.join(outpath, checksumf), 'w') as f2:
        for item in checksums:
            print('%-60s\t%s' % item, file=f2)


if __name__ == '__main__':

    #elvis = '6.5.X'
    # camptype = 'Full'  # Mini/Full
    #outpath = 'MiniCal_data02/Full_Camp_scripts_26JAN18_E6.5.0'
    # equipment = dict(operator='cpf',
    #                 sn_ccd1='D01',
    #                 sn_ccd2='D02',
    #                 sn_ccd3='D03',
    #                 sn_roe='SIM',
    #                 sn_rpsu='PSUs')
    # toWrite_def = OrderedDict(BIAS01=0, DARK01=0, CHINJ00=0, CHINJ01=0, CHINJ02=0,
    #                          FLAT01=0, FLAT02=0, PTC01=0, PTC02WAVE=0, PTC02TEMP=0, NL01=0,
    #                          PSF01=0, PSF02=0,
    #                          TP00=0, TP01=0, TP02=0,
    #                          PERSIST01=0, FOCUS00=0)

    # inputs = OrderedDict(camptype=camptype,elvis=elvis,outpath=outpath,
    #                     equipment=equipment,toWrite=toWrite)

    parser = OptionParser()
    parser.add_option("-j", "--json", dest="json",
                      default='', help="json file with inputs")

    (options, args) = parser.parse_args()

    if options.json == '':
        parser.print_help()
        sys.exit()

    inputs = vjson.load_jsonfile(options.json, useyaml=True)

    # MISSING: INPUTS VALIDATION

    if inputs['camptype'] == 'Full':
        test_generator = campaign.generate_test_sequence
    elif inputs['camptype'] == 'Mini':
        test_generator = minicampaign.generate_reduced_test_sequence
    elif inputs['camptype'] == 'MOT':
        test_generator = mot.generate_mot_sequence
    elif inputs['camptype'] == "devel":
        test_generator = devel.generate_test_sequence

    scwriter(inputs['toWrite'], test_generator, inputs['outpath'],
             inputs['equipment'], inputs['elvis'], inputs['CHAMBER'])
