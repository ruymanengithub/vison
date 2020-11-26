
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
from vison.support import cosmicrays
# END IMPORT

class EmulFPA(MetaPano):

    def __init__(self, testnames, **kwargs):

        super(EmulFPA, self).__init__(**kwargs)

        self.testnames = testnames

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

    def _f_tweak_quadrant(self, Qdata):
        """ """
        Qdata[-9:,:] = Qdata[-18:-9,:].copy() # filling up extended serial overscan
        return Qdata

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
                    
                    gain = self.cdps['GAIN'][block][CCDk][Q][0]

                    Qdata = kccdobj.get_quad(Q, canonical=True, extension=-1)

                    Qdata = self._f_tweak_quadrant(Qdata)

                    if CRs > 0:

                        if Qcounter == 0:
                            crImage = np.zeros_like(Qdata)
                            cosmics = cosmicrays.Cosmicrays(None, crImage, 
                                information=dict(exptime=CRexptime))
                        
                        Q_cr = cosmics.addToFluxTime(CRs, limit=None, verbose=False) / gain # ADU

                        Qdata[kccdobj.prescan:-kccdobj.overscan,vstart-1:vend] +=\
                            Q_cr[kccdobj.prescan:-kccdobj.overscan,vstart-1:vend]


                        Qdata[np.where(Qdata > 2**16-1)] = 2**16-1

                        cosmics.image = np.zeros_like(Qdata) # reset

                    kccdobj.set_quad(Qdata, Q,canonical=True, extension=-1)
                    Qcounter += 1


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

def run_TP11emul(CRs=0.,doLoad=False, doParse=False):
    """ """

    Tinputs = comminputs.copy()
    Tinputs.update(dict(
        testkey = 'TPX1',
        jsonf='ALLTESTS.json',
        respathroot=inpathroot,
        vcalfile = os.path.join(vcalpath,'CCD_CAL_CONVERSIONS_ALL_BLOCKS.xlsx')))


    emulator = EmulFPA(['TP01','TP11'], **Tinputs)
    emulator.prerun(doLoad=doLoad, doParse=doParse)

    stop() # how to deal with a test that changed named halfway through campaign?

    if CRs>0.:
        outfile = 'TPX1_emul_VGCCasFPA_CRs%.1f.fits' % CRs
    else:
        outfile = 'TPX1_emul_VGCCasFPA.fits'

    simulkwargs = dict(CRs=CRs,
        rep=1,
        relObsid=5,
        outfile=outfile,
        CRexptime=ROtime/2.)

    emulator.produce_emulation(**simulkwargs)

def run_TP21emul(CRs=0.,doLoad=False, doParse=False):
    """ """

def run_CHINJ01emul(CRs=0.,doLoad=False, doParse=False):
    """ """

    Cinputs = comminputs.copy()
    Cinputs.update(dict(
        testkey = 'CHINJ01',
        jsonf='ALLTESTS.json',
        respathroot=inpathroot,
        vcalfile = os.path.join(vcalpath,'CCD_CAL_CONVERSIONS_ALL_BLOCKS.xlsx')))


    emulator = EmulFPA(['CHINJ01'], **Cinputs)
    emulator.prerun(doLoad=doLoad, doParse=doParse)

    if CRs>0.:
        outfile = 'CHINJ01_emul_VGCCasFPA_CRs%.1f.fits' % CRs
    else:
        outfile = 'CHNJ01_emul_VGCCasFPA.fits'

    simulkwargs = dict(CRs=CRs,
        rep=1,
        relObsid=5,
        outfile=outfile,
        CRexptime=ROtime)

    emulator.produce_emulation(**simulkwargs)

def run_BIAS02emul(CRs=0.,doLoad=False, doParse=False):
    """ """

    Binputs = comminputs.copy()
    Binputs.update(dict(
        testkey = 'BIAS02',
        jsonf='ALLTESTS.json',
        respathroot=inpathroot,
        vcalfile = os.path.join(vcalpath,'CCD_CAL_CONVERSIONS_ALL_BLOCKS.xlsx')))


    emulator = EmulFPA(['BIAS02'], **Binputs)
    emulator.prerun(doLoad=doLoad, doParse=doParse)

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



def run_FFemul(CRs=0.,doLoad=False, doParse=False):
    """ """

    Finputs = comminputs.copy()
    Finputs.update(dict(
        testkey = 'FLAT02_730',
        jsonf='ALLTESTS.json',
        respathroot=inpathroot,
        vcalfile = os.path.join(vcalpath,'CCD_CAL_CONVERSIONS_ALL_BLOCKS.xlsx')))


    emulator = EmulFPA(['FLAT02_730'], **Finputs)
    emulator.prerun(doLoad=doLoad, doParse=doParse)

    if CRs>0.:
        outfile = 'FLAT02_730_emul_VGCCasFPA_CRs%.1f.fits' % CRs
    else:
        outfile = 'FLAT02_730_emul_VGCCasFPA.fits'

    simulkwargs = dict(CRs=CRs,
        rep=1,
        relObsid=5,
        outfile=outfile,
        CRexptime=ROtime/2.+10.)

    emulator.produce_emulation(**simulkwargs)




if __name__ == '__main__':

    doLoad = True
    doParse = True

    #run_CHINJ01emul(CRs=0.,doLoad=doLoad, doParse=doParse)

    #run_TP11emul(CRs=0.,doLoad=doLoad, doParse=doParse)
    #run_TP21emul(CRs=0.,doLoad=doLoad, doParse=doParse)

    #run_BIAS02emul(CRs=0.,doLoad=doLoad, doParse=doParse)
    #run_BIAS02emul(CRs=2.6,doLoad=doLoad, doParse=doParse)

    run_FFemul(CRs=0.,doLoad=doLoad, doParse=doParse)
    #run_FFemul(CRs=2.6,doLoad=doLoad, doParse=doParse)
