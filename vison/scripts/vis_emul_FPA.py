
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

    def tweak_quadrants(self, EMU, CRs=0., CRexptime=0., vstart=1,vend=2086):
        """
        Fixes additional serial overscan
        adds cosmic rays, if requested
        """

        Qcounter = 0

        for jY in range(1, EMU.fpaobj.fpamodel.NSLICES + 1):
            for iX in range(1, EMU.fpaobj.fpamodel.NSLICES + 1):
                CCDID = 'C_%i%i' % (jY, iX)
                
                locator = self.fpa.FPA_MAP[CCDID]

                block = locator[0]
                CCDk = locator[1]

                kccdobj = EMU.fpaobj.get_ccdobj(CCDID)

                for Q in self.Quads:
                    
                    gain = self.cdps['GAIN'][block][CCDk][Q]

                    Qdata = kccdobj.get_quad(Q, canonical=True, extension=-1)


                    Qdata[-9:,:] = Qdata[-18:-9,:].copy() # filling up extended serial overscan

                    if CRs > 0:

                        if ccounter == 0:
                            crImage = np.zeros_like(Qdata)
                            cosmics = cosmicrays(None, crImage, 
                                information=dict(exptime=CRexptime))

                        Q_cr = cosmics.addToFluxTime(CRs, limit=None, verbose=False) / gain

                        Qdata[kccdobj.prescan:kccdobj.overscan,vstart-1:vend] +=\
                            Q_cr[kccdobj.prescan:kccdobj.overscan,vstart-1:vend]

                        Qdata[np.where(Qdata > 2**16-1)] = 2**16-1

                        cosmics.image = np.zeros_like(Qdata) # reset


                    ccounter += 1

                    kccdobj.set_quad(Qdata, Q,canonical=True, extension=-1)


                EMU.fpaobj.set_ccdobj(kccdobj, CCDID, inextension=-1)

        return EMU




    def produce_emulation(self, **kwargs):
        """ """

        testname = self.testnames[0]
        irep = kwargs['rep']-1
        iobs = kwargs['relObsid']
        outputname = kwargs['outfile']
        CRs = kwargs['CRs']
        CRexptime = kwargs['CRexptime']

        ccd_dict = self.get_ccd_dict(testname,irep,iobs)

        EMUheader = dict(COSMICRA=CRs)

        EMU = cdpmod.LE1_CDP()

        EMU.ingest_inputs(ccd_dict, header=EMUheader, inextension=-1,
            fillval=0)

        EMU = self.tweak_quadrants(EMU, CRs, CRexptime)

        EMUname = os.path.join(self.outpathroot, outputname)

        EMU.savehardcopy(EMUname, clobber=True, uint16=True)
        print('Emulation saved in: %s' % EMUname)






outparent = 'FPA_SIMULS'
FPAdesign='final'
vcalpath = 'VOLT_CALIBS/ROE_VOLT_CALIBS'
inpathroot = 'FLIGHT'
comminputs = dict(outparent=outparent,
    design=FPAdesign,
    cdps=dict(gain=os.path.join('FPA_FINAL', 'PTC_FPA','GAIN_MX_PTC0X.pick')))
ROtime = 72. # seconds

def run_TP11emul():
    """ """

def run_TP21emul():
    """ """

def run_BIAS02emul(CRs=0.):
    """ """

    Binputs = comminputs.copy()
    Binputs.update(dict(
        testkey = 'BIAS02',
        jsonf='ALLTESTS.json',
        inpathroot=inpathroot,
        vcalfile = os.path.join(vcalpath,'CCD_CAL_CONVERSIONS_ALL_BLOCKS.xlsx')))


    emulator = EmulFPA('BIAS02', **Binputs)
    emulator.prerun(doLoad=False, doParse=False)

    if CRs>0.:
        outfile = 'BIAS02_emul_VGCCasFPA_CRs%.1f.fits' % CRs
    else:
        outfile = 'BIAS02_emul_VGCCasFPA.fits'

    simulkwargs = dict(CRs=CRs,
        rep=1,
        relObsid=5,
        outfile=outfile,
        CRexptime=ROtime/2.)

    emulator.produce_emulation(**simulkwargs)





def run_FFemul():
    """ """



if __name__ == '__main__':



    run_BIAS02emul(CRs=2.6)

