#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Classes and functions to generate ELVIS commanding scripts automatically.

Created on Wed May 24 15:31:54 2017

:author: Ruyman Azzollini
:contact: r.azzollini__at__ucl.ac.uk

"""

# IMPORT STUFF
from pdb import set_trace as stop
import pandas as pd
import string
import numpy as np
from collections import OrderedDict
import copy

from vison.support import context
#from vison.pipe import lib as pilib
# END IMPORT


script_dictionary = {'6.0.0': {
    'keys': ['frames',
             'program', 'test', 'IDL', 'IDH',
             'IG1_1_T', 'IG1_2_T', 'IG1_3_T',
             'IG1_1_B', 'IG1_2_B', 'IG1_3_B',
             'IG2_T', 'IG2_B',
             'OD_1_T', 'OD_1_B', 'OD_2_T', 'OD_2_B', 'OD_3_T', 'OD_3_B',
             'RD_T', 'RD_B',
             'iphi1', 'iphi2', 'iphi3', 'iphi4',
             'readmode_1', 'readmode_2',
             'ser_tpump', 'ser_tpump_mode',
             'vertical_clk', 'serial_clk',
             'flushes', 'exptime', 'shutter', 'electroshutter', 'vstart', 'vend',
             'sinvflush', 'sinvflushp',
             'chinj', 'chinj_rows_on', 'chinj_rows_off',  # 'chinj_repeat',
             'id_width', 'id_delay', 'chinj_ser_wait',
             'tpump', 'ser_shuffles', 'ver_shuffles',
             'dwell_v', 'dwell_h', 'tpump_mode',
             'motor', 'matrix_size', 'step_size',
             'add_h_overscan', 'add_v_overscan',
             'toi_flush', 'toi_tpump', 'toi_rdout', 'toi_chinj',
             'wavelength', 'pos_cal_mirror',
             'operator',
             'sn_ccd1', 'sn_ccd2', 'sn_ccd3', 'sn_roe', 'sn_rpsu',
             'comments'],

    # 'alias':dict(frames='Frames',
    # program='Program',test='Test',IDL='IDL(mV)',IDH='IDH(mV)',
    # IG1_1_T='IG1_1_T(mV)',IG1_2_T='IG1_2_T(mV)',IG1_3_T='IG1_3_T(mV)',
    # IG1_1_B='IG1_1_B(mV)',IG1_2_B='IG1_2_B(mV)',IG1_3_B='IG1_3_B(mV)',
    # IG2_T='IG2_T(mV)',IG2_B='IG2_B(mV)',
    # OD_1_T='OD CCD1 T(mV)',OD_1_B='OD CCD1 B(mV)',OD_2_T='OD CCD2 T(mV)',
    # OD_2_B='OD CCD2 B(mV)',OD_3_T='OD CCD3 T(mV)',OD_3_B='OD CCD3 B(mV)',
    # RD_T='RD T(mV)',RD_B='RD B(mV)',
    # iphi1='Iphi1',iphi2='Iphi2',iphi3='Iphi3',iphi4='Iphi4',
    # readmode_1='Readout mode R01',readmode_2='Readout mode R02',
    # ser_tpump='Serial T Pump Rdout',ser_tpump_mode='Serial T Pump mode',
    # vertical_clk='Vertical clk',serial_clk='Serial clk',
    # flushes='Flushes',exptime='Exposure time(ms)',
    # shutter='Shutter',electroshutter='Electronic shutter',
    # vstart='Vstart',vend='Vend',
    # sinvflush='S. Inv.Flush',sinvflushp='S.Inv.FlushPeriod',
    # chinj='charge injection',chinj_rows_on='chrg_inj rows on',
    # chinj_rows_off='chrg_inj rows off',#chinj_repeat='chrg_inj repeat',
    # id_width='id pulse length(us)',id_delay='id pulse delay(us)',
    # chinj_ser_wait='chinj_ser_wait',
    # tpump='Trap pumping',ser_shuffles='Serial shuffles',ver_shuffles='Vertical shuffles',
    # dwell_v='Dwell_V',dwell_h='Dwell_H',tpump_mode='T.Pump mode',
    # motor='Move motor',matrix_size='Matrix size',step_size='Step size',
    # add_h_overscan='Additional H.Overscan',add_v_overscan='Additional V.Overscan',
    # toi_flush='TOI Flush(us)',toi_tpump='TOI T.Pump(us)',
    #                    toi_rdout='TOI Rdout(us)',toi_chinj='TOI Chinj(us)',
    # wavelength='Wavelength',pos_cal_mirror='Pos cal mirror(mm)',
    #operator='Operator name',
    # sn_ccd1='CCD1 type/sn',sn_ccd2='CCD2 type/sn',sn_ccd3='CCD3 type/sn',
    # sn_roe='ROE type/sn',sn_rpsu='RPSU type/sn',
    # comments='Comments'),

    'defaults': dict(frames=1,
                     program='CALCAMP', test='Test', IDL=13000, IDH=18000,
                     IG1_1_T=6000, IG1_2_T=6000, IG1_3_T=6000,
                     IG1_1_B=6000, IG1_2_B=6000, IG1_3_B=6000,
                     IG2_T=5000, IG2_B=5000,
                     OD_1_T=26000, OD_1_B=26000, OD_2_T=26000, OD_2_B=26000, OD_3_T=26000, OD_3_B=26000,
                     RD_T=16000, RD_B=16000,
                     iphi1=1, iphi2=1, iphi3=1, iphi4=0,
                     readmode_1='Normal', readmode_2='Normal',
                     ser_tpump=0, ser_tpump_mode='S1&2',
                     vertical_clk='?', serial_clk='?',
                     flushes=7, exptime=0,
                     shutter='Thorlabs SC10', electroshutter=0,
                     vstart=0, vend=2066,
                     sinvflush=0, sinvflushp=500,
                     chinj=0, chinj_rows_on=1, chinj_rows_off=1,  # chinj_repeat=1,
                     id_width=1, id_delay=1, chinj_ser_wait=1,
                     tpump=0, ser_shuffles=0, ver_shuffles=0,
                     dwell_v=0, dwell_h=0, tpump_mode='V1&2',
                     motor=0, matrix_size=2, step_size=100,
                     add_h_overscan=0, add_v_overscan=0,
                     toi_flush=143, toi_tpump=1000, toi_rdout=1000, toi_chinj=1000,
                     wavelength='Filter 4', pos_cal_mirror=70,
                     operator='cpf',
                     sn_ccd1='CCD273-XX-X-XXX', sn_ccd2='CCD273-XX-X-XXX', sn_ccd3='CCD273-XX-X-XXX',
                     sn_roe='ROEXX', sn_rpsu='RPSUXX',
                     comments='')
}}

script_dictionary['6.1.0'] = script_dictionary['6.0.0']
script_dictionary['6.1.0']['keys'] = ['frames',
                                      'program', 'test', 'IDL', 'IDH',
                                      'IG1_1_T', 'IG1_2_T', 'IG1_3_T',
                                      'IG1_1_B', 'IG1_2_B', 'IG1_3_B',
                                      'IG2_T', 'IG2_B',
                                      'OD_1_T', 'OD_1_B', 'OD_2_T', 'OD_2_B', 'OD_3_T', 'OD_3_B',
                                      'RD_T', 'RD_B',
                                      'iphi1', 'iphi2', 'iphi3', 'iphi4',
                                      'readmode_1', 'readmode_2',
                                      'ser_tpump', 'ser_tpump_mode',
                                      'vertical_clk', 'serial_clk',
                                      'flushes', 'exptime', 'shutter', 'electroshutter', 'vstart', 'vend',
                                      'sumwell', 'inisweep',
                                      'sinvflush', 'sinvflushp',
                                      'chinj', 'chinj_rows_on', 'chinj_rows_off',  # 'chinj_repeat',
                                      'id_width', 'id_delay', 'chinj_ser_wait',
                                      'tpump', 'ser_shuffles', 'ver_shuffles',
                                      'dwell_v', 'dwell_h', 'tpump_mode',
                                      'motor', 'matrix_size', 'step_size',
                                      'add_h_overscan', 'add_v_overscan',
                                      'toi_flush', 'toi_tpump', 'toi_rdout', 'toi_chinj',
                                      'wavelength', 'pos_cal_mirror',
                                      'operator',
                                      'sn_ccd1', 'sn_ccd2', 'sn_ccd3', 'sn_roe', 'sn_rpsu',
                                      'comments']
# script_dictionary['6.1.0']['alias'].update(dict(sumwell='sumwell',
#                 inisweep='inisweep'))
script_dictionary['6.1.0']['defaults'].update(dict(sumwell=0,
                                                   inisweep=1))

script_dictionary['6.3.0'] = script_dictionary['6.1.0']
script_dictionary['6.3.0']['keys'] = [
    'frames', 'program', 'test', 'IDL', 'IDH',
    'IG1_1_T', 'IG1_2_T', 'IG1_3_T', 'IG1_1_B', 'IG1_2_B', 'IG1_3_B', 'IG2_T',
    'IG2_B', 'OD_1_T', 'OD_1_B', 'OD_2_T', 'OD_2_B', 'OD_3_T', 'OD_3_B',
    'RD_T', 'RD_B', 'IPHI1', 'IPHI2', 'IPHI3', 'IPHI4',
    'rdmode', 'flushes', 'siflsh', 'siflsh_p', 'swellw', 'swelldly', 'inisweep',
    'vstart', 'vend', 'toi_fl', 'toi_tp', 'toi_ro', 'toi_ch',
    'chinj', 'chinj_on', 'chinj_of', 'id_wid', 'id_dly', 'chin_dly',
    's_tpump', 's_tpmod', 'v_tpump', 'v_tpmod', 's_tp_cnt', 'v_tp_cnt', 'dwell_v',
    'dwell_s', 'exptime', 'shuttr', 'e_shuttr', 'mirr_on', 'mirr_pos', 'wave',
    'motr_on', 'motr_cnt', 'motr_siz', 'source',
    'operator', 'sn_ccd1', 'sn_ccd2', 'sn_ccd3', 'sn_roe', 'sn_rpsu',
    'comments']

script_dictionary['6.3.0']['defaults'] = {
    'frames': 1, 'program': 'CALCAMP', 'test': 'Test', 'IDL': 13.0, 'IDH': 18.0,
    'IG1_1_T': 6.0, 'IG1_2_T': 6.0, 'IG1_3_T': 6.0, 'IG1_1_B': 6.0,
    'IG1_2_B': 6.0, 'IG1_3_B': 6.0, 'IG2_T': 5.0, 'IG2_B': 5.0,
    'OD_1_T': 26.0, 'OD_1_B': 26.0, 'OD_2_T': 26.0, 'OD_2_B': 26.0,
    'OD_3_T': 26.0, 'OD_3_B': 26.0,
    'RD_T': 16.0, 'RD_B': 16.0,
    'IPHI1': 1, 'IPHI2': 1, 'IPHI3': 1, 'IPHI4': 0,
    'rdmode': 'fwd_bas',
    'flushes': 7, 'siflsh': 0, 'siflsh_p': 500,
    'swellw': context.sumwell['fwd_bas'][0], 'swelldly': context.sumwell['fwd_bas'][1],
    'inisweep': 1,
    'vstart': 0, 'vend': 2086,
    'toi_fl': 143, 'toi_tp': 1000, 'toi_ro': 1000, 'toi_ch': 1000,
    'chinj': 0, 'chinj_on': 1, 'chinj_of': 1,
    'id_wid': 100, 'id_dly': 1500, 'chin_dly': 0,
    's_tpump': 0, 's_tpmod': 23,
    'v_tpump': 0, 'v_tpmod': 123, 's_tp_cnt': 1000, 'v_tp_cnt': 1000,
    'dwell_v': 0, 'dwell_s': 0,
    'exptime': 0, 'shuttr': 1,
    'e_shuttr': 0,
    'mirr_on': 0, 'mirr_pos': 70,
    'wave': 4,
    'motr_on': 0, 'motr_cnt': 2, 'motr_siz': 50, 'source': 'flat',
    'operator': 'cpf',
    'sn_ccd1': 'CCD273-XX-X-XXX', 'sn_ccd2': 'CCD273-XX-X-XXX',
    'sn_ccd3': 'CCD273-XX-X-XXX', 'sn_roe': 'ROEXX', 'sn_rpsu': 'RPSUXX',
    'comments': ''
}

script_dictionary['6.5.X'] = dict()

script_dictionary['6.5.X']['keys'] = [
    'frames', 'program', 'test', 'IDL', 'IDH',
    'IG1_1_T', 'IG1_2_T', 'IG1_3_T', 'IG1_1_B', 'IG1_2_B', 'IG1_3_B', 'IG2_T',
    'IG2_B', 'OD_1_T', 'OD_1_B', 'OD_2_T', 'OD_2_B', 'OD_3_T', 'OD_3_B',
    'RD_T', 'RD_B', 'IPHI1', 'IPHI2', 'IPHI3', 'IPHI4',
    'rdmode', 'flushes', 'siflsh', 'siflsh_p', 'swellw', 'swelldly', 'inisweep',
    'vstart', 'vend', 'toi_fl', 'toi_tp', 'toi_ro', 'toi_ch',
    'chinj', 'chinj_on', 'chinj_of', 'id_wid', 'id_dly', 'chin_dly',
    's_tpump', 's_tpmod', 'v_tpump', 'v_tpmod', 's_tp_cnt', 'v_tp_cnt', 'dwell_v',
    'dwell_s', 'exptime', 'shuttr', 'e_shuttr', 'mirr_on', 'mirr_pos', 'wave',
    'motr_on', 'motr_cnt', 'motr_siz', 'source',
    'operator', 'sn_ccd1', 'sn_ccd2', 'sn_ccd3', 'sn_roe', 'sn_rpsu',
    'comments', 'roe-tab']
script_dictionary['6.5.X']['defaults'] = {
    'frames': 1, 'program': 'CALCAMP', 'test': 'Test',
    'IDL': 13.0, 'IDH': 18.0,
    'IG1_1_T': 6.0, 'IG1_2_T': 6.0, 'IG1_3_T': 6.0,
    'IG1_1_B': 6.0, 'IG1_2_B': 6.0, 'IG1_3_B': 6.0,
    'IG2_T': 5.0, 'IG2_B': 5.0,
    'OD_1_T': 26.0, 'OD_1_B': 26.0,
    'OD_2_T': 26.0, 'OD_2_B': 26.0,
    'OD_3_T': 26.0, 'OD_3_B': 26.0,
    'RD_T': 16.0, 'RD_B': 16.0,
    'IPHI1': 1, 'IPHI2': 1, 'IPHI3': 1, 'IPHI4': 0,
    'rdmode': 'fwd_bas',
    'flushes': 7, 'siflsh': 1, 'siflsh_p': 500,
    'swellw': context.sumwell['fwd_bas'][0], 'swelldly': context.sumwell['fwd_bas'][1],
    'inisweep': 1,
    'vstart': 0, 'vend': 2086,
    'toi_fl': 143, 'toi_tp': 1000, 'toi_ro': 1000, 'toi_ch': 1000,
    'chinj': 0, 'chinj_on': 1, 'chinj_of': 1,
    'id_wid': 100, 'id_dly': 1500, 'chin_dly': 0,
    's_tpump': 0, 's_tpmod': 23,
    'v_tpump': 0, 'v_tpmod': 123,
    's_tp_cnt': 1000, 'v_tp_cnt': 1000,
    'dwell_v': 0, 'dwell_s': 0,
    'exptime': 0, 'shuttr': 1, 'e_shuttr': 0,
    'mirr_on': 0, 'mirr_pos': 70,
    'wave': 4,
    'motr_on': 0, 'motr_cnt': 2, 'motr_siz': 50,
    'source': 'flat',
    'operator': 'cpf',
    'sn_ccd1': 'CCD273-XX-X-XXX', 'sn_ccd2': 'CCD273-XX-X-XXX',
    'sn_ccd3': 'CCD273-XX-X-XXX', 'sn_roe': 'ROEXX', 'sn_rpsu': 'RPSUXX',
    'comments': '',
    'roe-tab': 'hrdstart'
}

script_dictionary['7.2.X'] = copy.deepcopy(script_dictionary['6.5.X'])


def update_structdict(sdict, commvalues, diffvalues):
    """Updates an script structure with common values and differential values.

    :param sdict: dict, dictionary with script structure. Takes precedence over
                  commvalues.
    :param commvalues: dict, dictionary with common values to update sdict.
    :param diffvalues: dict, dictionaty with "differential" values to update
                       "sdict". Takes precedence over sdict and commvalues.

    """

    Ncols = sdict['Ncols']

    Ndiff = len(diffvalues.keys())

    for ic in range(1, Ncols+1):
        ickey = 'col%i' % ic

        for comkey in commvalues.keys():

            if comkey not in sdict[ickey].keys():
                sdict[ickey][comkey] = commvalues[comkey]

        if Ndiff > 0:
            sdict[ickey].update(diffvalues)

    return sdict


class Script(object):
    """Core Class that provides automatic test script generation and validation."""

    def __init__(self, defaults={}, structure={}, elvis=context.elvis):
        """Initialization.

        :param defaults: dict, default values.
        :param structure: dict, structure of script.
        :param elvis: ELVIS version.

        """

        self.scriptname = ''
        self.defaults = defaults
        self.structure = structure
        self.elvis = elvis
        self.cargo = []

    def build_cargo(self):
        """Updates 'cargo' attribute.
        'cargo': list of lists, each corresponding to a column in the script.
                 Each element in the inner lists is a register value.
                 The first column corresponds to the column with key names.

        Note: the number of frames is accumuled across columns, as ELVIS expects.

        """
        dictio = script_dictionary[self.elvis]
        keys = dictio['keys']
        #aliases = dictio['alias']

        stru = self.structure
        defaults = self.defaults

        self.cargo = []

        # self.cargo.append([aliases[key] for key in keys]) # first column, aliases
        self.cargo.append(keys)  # first column, aliases

        cols = ['col%i' % (i+1,) for i in range(stru['Ncols'])]

        cumframes = 0

        for col in cols:

            vallist = []

            for key in keys:
                if key in stru[col]:
                    val = stru[col][key]
                else:
                    val = defaults[key]
                vallist.append(val)

            frs = vallist[0]
            vallist[0] += cumframes  # accumulation of frames number
            cumframes += frs

            self.cargo.append(vallist)

        if self.elvis >= '6.0.0':
            endcol = ['']*len(keys)
            endcol[0] = 'End'
            self.cargo.append(endcol)

    def write(self, scriptname):
        """Writes self.cargo (script) to an excel file.

        :param scriptname: char, name of file where to write script.

        """

        ncols = len(self.cargo)

        colnames = ['col%04i' % i for i in range(ncols)]

        datadict = {}

        for ixc in range(ncols):
            datadict[colnames[ixc]] = self.cargo[ixc]

        df = pd.DataFrame(datadict)

        writer = pd.ExcelWriter(scriptname, engine='xlsxwriter')
        df.to_excel(writer, sheet_name='Sheet1', header=False, index=False)

        writer.save()

        self.scriptname = scriptname

        return None

    def get_struct_from_cargo(self):
        """ """
        Ncols = len(self.cargo[1:])

        keys = self.cargo[0]

        structure = dict()
        structure['Ncols'] = Ncols

        for i in range(1, Ncols+1):
            structure['col%i' % i] = dict(zip(keys, self.cargo[i]))

        return structure

    def load(self, *args, **kwargs):
        """alias method. Points to 'load_to_cargo'."""
        self.load_to_cargo(*args, **kwargs)

    def load_to_cargo(self, scriptname, elvis=context.elvis):
        """Loads an script from an excel file.

        :param scriptname: char, script to load
        :param elvis: char, ELVIS version of script to load

        """
        self.elvis = elvis

        converters = {'frames': str}
        df = pd.read_excel(scriptname, header=0, converters=converters)

        cargo = []
        cargo.append(['frames']+[string.strip(item)
                                 for item in df['frames'].tolist()])

        sdict = script_dictionary[self.elvis]
        keys = sdict['keys']
        #alias = sdict['alias']
        #aliaslist = [alias[key] for key in keys]

        assert np.all(cargo[0] == keys)

        for ixc in df.columns:
            if ixc == 'frames':
                continue
            elif ixc == 'End':
                break
            elif isinstance(ixc, int):
                rawvals = [ixc] + df[ixc].tolist()
                vals = []
                for item in rawvals:
                    if isinstance(item, unicode):
                        vals.append(string.strip(str(item)))
                    else:
                        vals.append(item)

                cargo.append(vals)

        self.cargo = cargo

        self.scriptname = scriptname

    def validate(self, defaults, structure, elvis=context.elvis):
        """Not sure 'validation' will work like as implemented...
        TODO: validate self.validate"""

        assert self.elvis == elvis

        access = Script(defaults, structure, elvis=elvis)
        access.build_cargo()

        assert len(access.cargo) == len(self.cargo)

        for i in range(len(self.cargo)):
            assert np.all(self.cargo[i] == access.cargo[i])


def test0():
    """ """

    scriptname = '../../../../CALIBRATION/CAMPAIGN/ELVIS/ELVIS_6.3.0_27Jan17/EUCLID_QM_ELVIS_script_v6.3.0.xlsx'

    script = Script()
    script.load(scriptname)
    stop()


if __name__ == '__main__':

    test0()
