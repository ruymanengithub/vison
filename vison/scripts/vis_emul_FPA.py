
# IMPORT STUFF
from pdb import set_trace as stop
import os
import numpy as np
import copy

from vison.metatests.pano import MetaPano
from vison.metatests.bias import MetaBias
from vison.support import cosmicrays as crs
from vison.support import files
from vison.datamodel import ccd as ccdmod
from vison.datamodel import cdp as cdpmod
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


    def get_ccd_dict(self, testname, irep, iobs):
        """ """

        #PT = self.ParsedTable[testname]

        ccd_dict = dict()
        CCDs = ['CCD1','CCD2','CCD3']

        for jY in range(self.NCOLS_FPA):
            for iX in range(self.NSLICES_FPA):
                Ckey = 'C_%i%i' % (jY + 1, iX + 1)

                locator = self.fpa.FPA_MAP[Ckey]

                block = locator[0]
                CCDk = locator[1]
                jCCD = CCDs.index(CCDk)

                inventoryitem = self.inventory[block][testname][irep]

                iFileName = inventoryitem['dd'].mx['File_name'][iobs,jCCD]
                ireldatapath = inventoryitem['dd'].mx['datapath'][iobs,jCCD]

                idatapath = os.path.join('FLIGHT',block,\
                    *ireldatapath.split(os.path.sep)[1:])

                iFITS = os.path.join(idatapath,'{}.fits'.format(iFileName))

                iccdobj = ccdmod.CCD(iFITS)

                ccd_dict[Ckey] = copy.deepcopy(iccdobj)
        
        return ccd_dict


                


    def produce_emulation(self, **kwargs):
        """ """

        testname = self.testnames[0]
        irep = kwargs['rep']-1
        iobs = kwargs['relObsid']
        outputname = kwargs['outfile']
        

        ccd_dict = self.get_ccd_dict(testname,irep,iobs)

        EMUheader = dict()


        EMUcdp = cdpmod.LE1_CDP()

        EMUcdp.ingest_inputs(ccd_dict, header=EMUheader, inextension=-1,
            fillval=0)

        EMUcdpname = outputname

        EMUcdp.savehardcopy(EMUcdpname, clobber=True, uint16=True)
        print('Emulation saved in: %s' % EMUcdpname)






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
        relObsid=5,
        outfile=os.path.join(respathroot,'BIAS02_emul_VGCCasFPA.fits'))

    emulator.produce_emulation(**simulkwargs)





def run_FFemul():
    """ """



if __name__ == '__main__':



    run_BIAS02emul(addCRs=True)

