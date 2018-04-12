#! /usr/bin/env python

"""

Just a collection of LaTeX templates for use in report.py

:History:
Created on Mon Jan 30 2017

:author: Ruyman Azzollini 
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
from pdb import set_trace as stop
import os
from vison import data
import string as st
#import numpy as num
#from PIL import Image
#import string as st

# END IMPORT

isthere = os.path.exists


def generate_header(test, model, author):
    """ """

    headertexf = os.path.join(data.__path__[0], 'header_template.tex')

    fp = open(headertexf, 'r')
    headertexlines = fp.readlines()
    fp.close()

    headertexlines = [item[0:-1] for item in headertexlines]  # removes end \n

    header = st.join(headertexlines, 'JOIN___LINES')
    header = header.replace('%', '%%')
    header = header.replace('__PLACEHOLDER__', '%s')
    header = header % (author, model, test)
    header = header.replace('%%', '%')

    headerList = st.split(header, 'JOIN___LINES')

    return headerList


def generate_preamble(model, test, custodian='Ruyman Azzollini'):

    preambletexf = os.path.join(data.__path__[0], 'preamble_template.tex')

    fp = open(preambletexf, 'r')
    preambletexlines = fp.readlines()
    fp.close()

    preambletexlines = [item[0:-1]
                        for item in preambletexlines]  # removes end \n

    preamble = st.join(preambletexlines, 'JOIN___LINES')
    preamble = preamble.replace('%', '%%')
    preamble = preamble.replace('__PLACEHOLDER__', '%s')
    preamble = preamble % (model, test, custodian)
    preamble = preamble.replace('%%', '%')

    preambleList = st.split(preamble, 'JOIN___LINES')

    return preambleList
