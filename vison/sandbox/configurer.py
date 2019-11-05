#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 17:14:00 2018

:author: Ruyman Azzollini

"""


class Config():
    """ """

    def __init__(self, calpath, resroot, BLOCKID, inCDPs, equipment, elvis):
        """ """
        self.calpath = calpath
        self.resroot = resroot
        self.elvis = elvis
        self.BLOCKID = BLOCKID
        self.inCDPs = inCDPs
        self.equipment = equipment
        self.elvis = elvis

        self.inputsdict = dict()

    def init_inputsdict(self, tasks):
        """ """
        self.inputsdict = dict(BLOCKID=self.BLOCKID,
                               resultsroot=self.resroot,
                               tasks=[])

    def get_COMMTASK_inp(self, OBSID_lims, datapath, explogf):
        inputs = dict(
            OBSID_lims=OBSID_lims,
            datapath=datapath,
            explogf=explogf,
            elvis=self.elvis,
            diffvalues=self.equipment,
            inCDPs=self.inCDPs
        )
        return inputs

    def get_BIAS01_inp(self, N, OBSID_lims, datapath, explogf):
        """ """
        inputs = self.get_COMMTASK_inp(OBSID_lims, datapath, explogf)

        inputs.update(
            test='BIAS01',
            N=N,
            resultspath='BIAS01',
            todo_flags=dict(init=True, check=True, prep=True,
                            basic=True, meta=True, report=True),
        )
        return inputs

    # def set_BIAS01_inp(self,*args,**kwargs):
    #    inps = self.get_BIAS01_inp(self,*args,**kwargs)
    #    self.inputsdict['BIAS01'] = inps

    def get_DARK01_inp(self, OBSID_lims, datapath, explogf):
        inputs = self.get_COMMTASK_inp(OBSID_lims, datapath, explogf)
        inputs.update(
            test='DARK01',
            resultspath='DARK01',
            todo_flags=dict(init=True, check=True, prep=True, basic=True, meta=True, report=True))
        return inputs

    def get_CHINJ01_inp(self, OBSID_lims, datapath, explogf):
        inputs = self.get_COMMTASK_inp(OBSID_lims, datapath, explogf)
        toi_chinj01 = 500
        inputs.update(
            test='CHINJ01',
            resultspath='CHINJ01',
            todo_flags=dict(init=True, check=True, prep=True,
                            basic=True, meta=True, report=True),
            IDL=11.,
            IDH=18.,
            IG1s=[2., 6.],
            id_delays=[toi_chinj01 * 3, toi_chinj01 * 2],
            toi_chinj=toi_chinj01)
        return inputs

    def set_TASK_inp(self, *args, **kwargs):
        taskname = args[0]
        builder = self.taskdict[taskname]
        inps = builder(self, *args[1:], **kwargs)
        self.inputsdict[taskname] = inps
