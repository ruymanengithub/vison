
# IMPORT STUFF
from pdb import set_trace as stop
import os
import numpy as np

from vison.metatests.pano import MetaPano
from vison.metatests.bias import MetaBias
from vison.support import cosmicrays as crs
from vison.support import files
# END IMPORT

class EmulFPA(MetaPano):

    def __init__(self, testname, **kwargs):

        super(EmulFPA, self).__init__(**kwargs)

        self.testnames = [testname]

    def prerun(self, doLoad=True, doParse=True, doEmul=True):
        """ """

        if not os.path.exists(self.outpathroot):
            os.system('mkdir -p %s' % self.outpathroot)

        dictiopick = os.path.join(self.outpathroot,
                                  '%s_dictionary.pick' % self.testkey.lower())


        if doLoad:
            print('Loading block(s) results...')
            self.load_block_results()
            files.cPickleDumpDictionary(self.inventory, dictiopick)
        else:
            print('Re-loading block(s) results')
            self.inventory = files.cPickleRead(dictiopick)

        parsedpick = os.path.join(self.outpathroot,
                                  '%s_parsed.pick' % self.testkey.lower())

        if doParse:
            print('Parsing results')
            for testname in self.testnames:
                self.parse_test_results(testname)
            parsedbundle = dict(PT=self.ParsedTable,
                                products=self.products)
            files.cPickleDumpDictionary(parsedbundle, parsedpick)

        else:
            print('Re-loading parsed results')
            parsedbundle = files.cPickleRead(parsedpick)
            self.ParsedTable = parsedbundle['PT'].copy()
            self.products = parsedbundle['products'].copy()


    def produce_emulation(self, **kwargs):
        """ """

        stop()


outparent = 'FPA_SIMULS'
FPAdesign='final'
vcalpath = 'VOLT_CALIBS/ROE_VOLT_CALIBS'
respathroot = 'FLIGHT'
comminputs = dict(outparent=outparent,
    design=FPAdesign,
    cdps=dict(gain=os.path.join('FPA_FINAL', 'PTC_FPA','GAIN_MX_PTC0X.pick')))


def run_TP11emul():
    """ """

def run_TP21emul():
    """ """

def run_BIAS02emul(addCRs=False):
    """ """

    Binputs = comminputs.copy()
    Binputs.update(dict(
        testkey = 'BIAS02',
        jsonf='ALLTESTS.json',
        respathroot=respathroot,
        vcalfile = os.path.join(vcalpath,'CCD_CAL_CONVERSIONS_ALL_BLOCKS.xlsx')))


    emulator = EmulFPA('BIAS02', **Binputs)
    emulator.prerun(doLoad=False, doParse=False)

    simulkwargs = dict(addCRs=addCRs,
        rep=1,
        relObsid=5)

    emulator.produce_emulation(**simulkwargs)





def run_FFemul():
    """ """



if __name__ == '__main__':



    run_BIAS02emul(addCRs=True)

