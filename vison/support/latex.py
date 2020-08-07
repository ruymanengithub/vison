#! /usr/bin/env python

"""

Just a collection of LaTeX-generating functions for use in report.py

:History:
Created on Mon Jan 30 2017

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import os
from vison import data
import string as st

# END IMPORT

isthere = os.path.exists


def replace_in_template(texf, values):
    """ """

    fp = open(texf, 'r')
    texlines = fp.readlines()
    fp.close()

    texlines = [item[0:-1] for item in texlines]  # removes end \n

    tex = st.join(texlines, 'JOIN___LINES')
    tex = tex.replace('%', '%%')
    tex = tex.replace('__PLACEHOLDER__', '%s')
    tex = tex % values
    tex = tex.replace('%%', '%')

    texList = st.split(tex, 'JOIN___LINES')

    return texList


def generate_header(test, model, author, reference='7-XXX', issue=0.0, doDraft=False):
    """ """

    if doDraft:
        headertexf = os.path.join(data.__path__[0], 'header_draft_template.tex')    
    else:
        headertexf = os.path.join(data.__path__[0], 'header_template.tex')

    headerList = replace_in_template(
        headertexf, (author, model, test, reference, issue))
    return headerList


def generate_preamble(model, test, custodian='Ruyman Azzollini', reference='7-XXX',
                      preambletexf=None):
    if preambletexf is None:
        preambletexf = os.path.join(data.__path__[0], 'preamble_template.tex')
    preambleList = replace_in_template(
        preambletexf, (model, test, reference, custodian))
    return preambleList
