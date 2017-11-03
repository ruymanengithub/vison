# -*- coding: utf-8 -*-
"""

Euclid-VIS Calibration Programme Pipeline: vison

Reporting Utilities.

:History:
Created on Wed Jan 25 16:58:33 2017

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
import sys
from pdb import set_trace as stop
import os
from vison.support import vistime
#from latex import generate_header, generate_preamble
import latex as lx
from vison import data as visondata

from astropy import table as astable
from astropy.io import ascii

import tempfile
import string as st
# END IMPORT




default_header = [
'\documentclass[12pt,a4paper]{article}',\
'\usepackage{t1enc}',\
'\usepackage{a4wide}',\
'\usepackage{latexsym}',\
'\usepackage[dvips]{graphicx}',\
'%\usepackage{psfig}',\
'\usepackage[spanish]{babel}',\
'\usepackage[ansinew]{inputenc}',\
'\usepackage{amsmath}\usepackage{amsfonts}\usepackage{amssymb}',\
'\usepackage{fancyhdr}',\
'\usepackage{longtable}',\
'\usepackage{multicol}',\
'\\addtolength{\marginparwidth}{-1cm}',\
'\\addtolength{\headwidth}{\marginparsep}',\
'\\addtolength{\headwidth}{\marginparwidth}',\
'\\addtolength{\\textwidth}{1.5cm}',\
'\n']

class Content(object):
    """ """
    
    def __init__(self,contenttype=''):
        """ """
        self.type=contenttype
    
    def __str__(self):
        """ """
        return self.type

class Container(Content,object):
    """ """
    
    def __init__(self):
        """ """
        self.Contents = []
        
    def __str__(self):
        """ """
        txt = ''
        for item in self.Contents:
            txt += self.Contents.__str__()
        return txt
    
    def add_to_Contents(self,item):
        """ """
        
        try:
            if len(self.Contents)>0:
                lastCont = self.Contents[-1]
                if not isinstance(lastCont,Section):
                    self.Contents.append(item)
                if isinstance(item,Section):
                    lastlevel = lastCont.level
                    itemlevel = item.level
                    if itemlevel > lastlevel:
                        self.Contents[-1].add_to_Contents(item)
                    else:
                        self.Contents.append(item)
                else:
                    self.Contents[-1].add_to_Contents(item)
            else:
                self.Contents.append(item)
        except AttributeError:
            self.Contents.append(item)
        
        return None

class Section(Container):
    """ """
    
    def __init__(self,keyword,Title='',level=0):
        
        self.levelnames = ['section','subsection','subsubsection']
        
        self.type = self.levelnames[level]
        self.level = level
        self.Title = Title
        self.keyword = keyword
        self.Contents = []
        
    def generate_Latex(self):
        """ """
        tex = ['\%s{%s}' % (self.type,self.Title),'\n']
        
        for item in self.Contents:            
            tex += item.generate_Latex()
        
        return tex
    
    def __str__(self):
        
        txt = '%s (%s) ' % (self.type, self.keyword)
        for item in self.Contents:
            txt += item.__str__()
        return txt

class Figure(Content):
    """ """
    def __init__(self,figpath,textfraction=0.7,caption=None,label=None):
        
        self.type = 'figure'
        self.figpath = figpath
        self.textfraction = textfraction
        self.caption = caption
        self.label = None
        
        
    def generate_Latex(self):
        """Generates LaTeX as list of strings."""
        
        
        tex_raw = ['\\begin{figure}[!htb]',
        '\centering',
        '\includegraphics[width=__PLACEHOLDER__\\textwidth]{__PLACEHOLDER__}',
        '\\vspace{-1mm}',
        '\end{figure}']
        
        substitutes = [self.textfraction,self.figpath]
        
        if self.caption is not None:
            tex_raw.insert(-1,'\caption{__PLACEHOLDER__}')
            substitutes.append(self.caption)
        
        if self.label is not None:
            tex_raw.insert(-1,'\label{__PLACEHOLDER__}')
            substitutes.append(self.label)
        
        texjoint = st.join(tex_raw,'JOIN_LINES')
        texjoint = texjoint.replace('__PLACEHOLDER__','%s')
        texjoint = texjoint % tuple(substitutes)
        
        tex = st.split(texjoint,'JOIN_LINES')
        tex += ['\n'] # just improves legigibility of .tex
        
        return tex


class Table(Content):
    """ 
    PENDING:
        - adjust width of table to texwidth:
            \resizebox{\textwidth}{!}{
            ...
            \end{tabular}}
        - include option to rotate table to show in landscape
    
    """
    def __init__(self,tableDict,formats=dict(),names=[],caption=None):
        
        table = astable.Table(data=tableDict,names=names)
        
        self.type = 'table'
        self.table = table
        self.formatDict = formats
        self.caption = caption
        
    def generate_Latex(self,):
        """Generates LaTeX as list of strings."""
        
        table = self.table
        formats = self.formatDict
        
        tf = tempfile.TemporaryFile()
        
        ascii.write(table,output=tf,formats=formats,Writer=ascii.Latex)
        
        tf.seek(0)
        tex = tf.readlines()
        tf.close()
        
        if self.caption is not None:
            captiontex = '\caption{%s}\n' % self.caption
            tex.insert(-1,captiontex)
        
        tex = [item[0:-1] for item in tex] # removes ending \n
        
        tex += ['\n'] # just improves legigibility of .tex
        
        
        tex = [st.replace(item,'_','\_') for item in tex]
        
        return tex
    
    
class Text(Content):
    """ """
    def __init__(self,text):
        
        self.type = 'text'
        self.text = text
        
    def generate_Latex(self):
        """ """
        if isinstance(self.text,str):
            tex = [self.text]
        elif isinstance(self.text,list):
            tex = self.text
        
        tex += ['\n'] # just improves legigibility of .tex
        
        return tex

class Report(Container):
    
    
    def __init__(self,TestName='Test',Model='XM',Contents=[],Texheader=default_header,
                 Texbody=['\\begin{document}']):
        """ """
        
        self.TestName = st.replace(TestName,'_','\_')
        self.Model = Model
        self.timestamp = vistime.get_time_tag()
        
        self.Texheader = Texheader
        self.Texbody = Texbody
        self.Texfooter = ['\end{document}']
        
        self.Contents = Contents
    
        
    
    def generate_Header(self,author='Ruyman Azzollini'):
        
        test = self.TestName
        model = self.Model
        
        headerList = lx.generate_header(test,model,author)        
        self.Texheader = headerList
        
    
    def generate_Texbody(self,custodian='Ruyman Azzollini'):
        """ """
        model = self.Model
        test = self.TestName
        
        Texbody = lx.generate_preamble(model,test,custodian)
        
        for item in self.Contents:            
            Texbody += item.generate_Latex()
                        
        self.Texbody = Texbody
    
    def has_Section(self,keyword):
        
        for item in self.Contents:
            if isinstance(item,Section):
                if item.keyword == keyword: return True
        
        return False
    
    def add_Section(self,keyword='',Title='',level=0,overwrite=True):
        """ """
        if overwrite and self.has_Section(keyword):
            self.drop_Section(keyword)
        self.add_to_Contents(Section(keyword,Title,level))
    
    def drop_Section(self,keyword):
        
        ix2remove = []
        
        Ncont = len(self.Contents)
        
        for ix,item in enumerate(self.Contents):
            if isinstance(item,Section):
                if item.keyword == keyword:
                    ix2remove.append(ix)
                    
                    children = range(ix,Ncont-1)
                    
                    if len(children) == 0: break
                    
                    for j in range(len(children)):
                        _item = self.Contents[children[j]]
                        if isinstance(_item,Section):
                            if _item.level > self.Contents[ix].level:
                                ix2remove.append(children[j])
                            else:
                                break
                        else:
                            break
                    
        for ix in ix2remove:
            self.Contents.pop(ix)
        
        
    
    def add_Figure(self,figpath,texfraction=0.7,caption=None,label=None):
        """ """
        self.add_to_Contents(Figure(figpath,texfraction,caption,label))
    
    
    def add_Table(self,tableDict,formats=dict(),names=[],caption=''):
        """ """
        # tableDict,formats=dict(),names=[],caption=None
        self.add_to_Contents(Table(tableDict,formats,names,caption))
    
    
    def add_text(self,text):
        """ """
        self.add_to_Contents(Text(text))
    
    
    def doreport(self,reportname,cleanafter=False):
        
        self.generate_Header()
        self.generate_Texbody()
        outfiles = self.writeto(reportname,cleanafter)
        return outfiles
        
    
    def writeto(self,fileroot,cleanafter=False):
        """Writes pdf with contents"""
        
        fileLaTex = '%s.tex' % fileroot
        
        if os.path.exists(fileLaTex):
            os.system('rm %s' % fileLaTex)
        
        self.writeLaTex(fileLaTex)
        outfiles = self.compileLaTex2pdf(fileroot,cleanafter=cleanafter)
        
        return outfiles
        
    def writeLaTex(self,filename):
        """Writes the LaTeX file"""
        
        f = open(filename,'a')
        for line in self.Texheader : print >> f, line
        for line in self.Texbody : print >> f, line
        for line in self.Texfooter: print >> f, line
        f.close()
        
    
    
    def compileLaTex2pdf(self,fileroot,cleanafter=False):
        """Compiles a Latex file"""
        
        EuclidViscls = 'EuclidVIS.cls'
        logo = 'logosmall.png'
        deluxetablesty = 'deluxetable.sty'
        
        os.system('ln -s %s' % os.path.join(visondata.__path__[0],EuclidViscls))
        os.system('ln -s %s' % os.path.join(visondata.__path__[0],logo))
        os.system('ln -s %s' % os.path.join(visondata.__path__[0],deluxetablesty))
        
        execline1 = 'latex %s.tex' % fileroot
        os.system(execline1)
        execline2 = 'dvipdf %s.dvi %s.pdf' % tuple([fileroot]*2)
        os.system(execline2)
        os.system(execline2) # twice to get all references
        
        
        if cleanafter :
            stop()
            os.system('rm %s.dvi %s.aux %s.log %s.tex %s.out %s.soc %s.toc' % \
            tuple([fileroot]*7))
            outfiles = [item % fileroot for item in \
              ['%s.pdf']]
        else:
            outfiles = [item % fileroot for item in \
              ['%s.pdf','%s.dvi','%s.aux','%s.log','%s.tex','%s.out','%s.soc',\
              '%s.toc']]
              
        os.system('rm %s' % EuclidViscls)
        os.system('rm %s' % logo)
        os.system('rm %s' % deluxetablesty)
        
        return outfiles
    
    
if __name__ == '__main__':
    
    report = Report(TestName='Dummy',Model='QM')
    
    outpath = './'
    if not os.path.exists(outpath):
        os.system('mkdir %s' % outpath)
    outfile = os.path.join(outpath,'Test_Report')
    
    report.add_Section('s1','Intro',level=0)
    report.add_text('this is a section')
    report.add_Section('s11','Sub-Intro',level=1)
    report.add_text('this is a sub-section')
    report.add_Section('s2','First Section',level=0)
    report.add_text('this is inside First Section')
    report.add_Section('s21','First Section Subsection',level=1)
    report.add_text('this is inside First Section Subsection')
    report.add_Section('s3','Second Section',level=0)
    report.add_text('this is inside Second Section')
    report.drop_Section('s2')
    
    report.doreport('Test_Report',cleanafter=False)