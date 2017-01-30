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
from pdb import set_trace as stop
import os
from vison.pipe.lib import get_time_tag
from latex import generate_header, generate_preamble
from vison import data as visondata
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
'\n',\
'\\begin{document}']


class Content(dict):
    """ """
    
    def __init__(self,contenttype=''):
        """ """
        self.type=contenttype
        

class Section(Content):
    """ """
    
    def __init__(self,Title='',level=0):
        
        levelnames = ['section','subsection','subsubsection']
        
        self.type = levelnames[level]
        self.Title = Title
        
    def generate_Latex(self):
        """ """
        tex = ['\%s{%s}' % (self.type,self.Title)]
        return tex

class Figure(Content):
    """ """
    def __init__(self,figpath,caption=''):
        
        self.type = 'figure'
        self.caption = caption
        
    def generate_Latex(self):
        """PENDING"""
        tex = ['']
        return tex


class Table(Content):
    """ """
    def __init__(self,tableDict,caption=''):
        
        self.type = 'table'
        self.caption = caption
        
    def generate_Latex(self):
        """PENDING"""
        tex = ['']
        return tex
    
    
class Text(Content):
    """ """
    def __init__(self,text):
        
        self.type = 'text'
        self.text = text
        
    def generate_Latex(self):
        """ """
        tex = self.text
        return tex

class Report(object):
    
    
    def __init__(self,TestName='Test',Model='XM',Contents=[],Texheader=default_header,Texbody=['']):
        """ """
        
        self.TestName = TestName
        self.Model = Model
        self.timestamp = get_time_tag()
        
        self.Texheader = Texheader
        self.Texbody = Texbody
        self.Texfooter = ['\end{document}']
        
        self.Contents = Contents
    
    
    def generate_Header(self,author='Ruyman Azzollini'):
        
        test = self.TestName
        model = self.Model
        
        headerList = generate_header(test,model,author)        
        self.Texheader = headerList
        
    
    def generate_Texbody(self,custodian='Ruyman Azzollini'):
        """ """
        model = self.Model
        test = self.TestName
        
        Texbody = generate_preamble(model,test,custodian)
        
        for item in self.Contents:            
            Texbody += item.generate_Texbody()
                        
        self.Texbody = Texbody
    
    def add_Section(self,Title='',level=0):
        """ """
        self.Contents.append(Section(Title,level))
        
        
    def add_Figure(self,figpath,caption=''):
        """ """
        self.Contents.append(Figure(figpath,caption))
        
    def add_Table(self,tableDict,caption=''):
        """ """
        self.Contents.append(Table(tableDict,caption))
        
    def add_text(self,text):
        """ """
        self.Contents.append(Text(text))
    
    
    def writeto(self,fileroot,cleanafter=False):
        """Writes pdf with contents"""
        
        fileLaTex = '%s.tex' % fileroot
        
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
        
        os.system('ln -s %s' % os.path.join(visondata.__path__[0],EuclidViscls))
        os.system('ln -s %s' % os.path.join(visondata.__path__[0],logo))
        
        execline1 = 'latex %s.tex' % fileroot
        os.system(execline1)
        execline2 = 'dvipdf %s.dvi %s.pdf' % tuple([fileroot]*2)
        os.system(execline2)
        os.system(execline2) # twice to get all references
        
        if cleanafter :
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

        return outfiles