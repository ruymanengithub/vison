#! /usr/bin/env python

"""Module with a class to create LaTeX files from python (scripting).

:History:
Created on Fri Mar 18 2016

:author: Ruyman Azzollini 
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
from pdb import set_trace as stop
import os
import numpy as num
from PIL import Image
import string as st

# END IMPORT

isthere = os.path.exists
latex_header = [
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
latex_footer = ['\end{document}']


class LaTeX(dict):
    """Class to write and compile a generic LaTeX file"""    
    
    def __init__(self,fontsize=12):
        self.body = []
        self.header = latex_header
        self.header[0] = self.header[0].replace('12pt','%ipt' % fontsize)        
        self.footer = latex_footer
    
    def DoBody(self,header,figures,figcomms={}):
        """Does a Body"""
        
        ImagesperRow = 3
        imgsize = 5
        
        ncols = int(num.ceil(len(figurelist)/float(ImagesperRow)) * 2)
        #ncells = ncols * ImagesperRow
        
        body = []
        body.append(header)
        body.append('\\begin{longtable}{|%s}' % ('c|' * ImagesperRow))
        body.append('\\hline')
        
        figcounter = 0
        commcounter = 0
        
        for j in range(ncols):
            for i in range(ImagesperRow):
                
                if j % 2 == 0:
                    try : idf = figurelist[figcounter]    
                    except IndexError: idf = 'Empty'
                    
                    try:
                        image = figures[idf]
                        size = Image.open(image).size
                        if size[0] >= size[1]:
                            sizekey = 'width'
                        elif size[1] > size[0]:
                            sizekey = 'height'
                        newrow = '\includegraphics[%s=%scm]{%s}&' % \
                        (sizekey,imgsize,image)
                    except KeyError:
                        idf = idf.replace('_','\_')
                        newrow = '$%s$   &' % idf
                    if (i+1) == ImagesperRow:
                        newrow = newrow[0:-1] + ' \\\\'
                    body.append(newrow)
                    if (i+1) % ImagesperRow ==0: body.append('\\hline')
                    figcounter +=1
                else : 
                    try: idc = figurelist[commcounter]    
                    except IndexError: idc = 'silent'
                    try:
                        comment = figcomms[idc]
                        newrow = '%s &' % (comment)
                    except KeyError:
                        newrow = '   &'
                    if (i+1) == ImagesperRow:
                        newrow = newrow[0:-1] + ' \\\\'
                    body.append(newrow)
                    if (i+1) % ImagesperRow ==0: body.append('\\hline')
                    commcounter +=1
        
        body.append('\end{longtable}')
        
        self.body = body
        return None
    
    def GenDoBody(self,header,figures,ImagesperRow=3,imgsize=5):
        """Does a Body"""        
        
                

        ncols = int(num.ceil(len(figures)/float(ImagesperRow)) * 2)
        #ncells = ncols * ImagesperRow
        
        body = []
        if header != '': body.append(header+'\\\\')
        body.append('\\begin{longtable}{|%s}' % ('c|' * ImagesperRow))
        body.append('\\hline')
        
        for j in range(ncols):
            for i in range(ImagesperRow):
                              
                try: image = figures[j*ImagesperRow+i]
                except IndexError: break
                size = Image.open(image).size
                if size[0] >= size[1]:
                    sizekey = 'width'
                elif size[1] > size[0]:
                    sizekey = 'height'
                newrow = '\includegraphics[%s=%scm]{%s}&' % \
                    (sizekey,imgsize,image)
    		     
                if (i+1) == ImagesperRow:
                    newrow = newrow[0:-1] + ' \\\\'
                body.append(newrow)
                if (i+1) % ImagesperRow ==0: body.append('\\hline')
                       
        body.append('\end{longtable}')
        
        self.body = body

        return None
    
    
    def addfigtobody(self,eps,imgsize=5):
        """ """
        
        size = Image.open(eps).size
        if size[0] >= size[1]:
            sizekey = 'width'
        elif size[1] > size[0]:
            sizekey = 'height'
        
        newline = '\includegraphics[%s=%scm]{%s}' % \
                    (sizekey,imgsize,eps)
        
        self.body.append(newline)
        
    
    def Write(self,filename):
        """Writes the LaTeX file"""
        
        f = open(filename,'a')
        for line in self.header : print >> f, line
        for line in self.body : print >> f, line
        for line in self.footer : print >> f, line
        f.close()
    
    def Compile(self,LatexFile,cleanafter=False,figures=[]):
        """Compiles a Latex file"""
        # IMPORT STUFF
        import os
        import string
        # END IMPORT
        
        root = LatexFile[0:string.rfind(LatexFile,'.tex')]
        
        execline1 = 'latex %s' % LatexFile
        os.system(execline1)
        execline2 = 'dvips -o %s %s' % (root+'.ps',root+'.dvi')
        os.system(execline2)
        
        if cleanafter :
            os.system('rm %s %s %s %s' % \
            (root+'.dvi',root+'.aux',root+'.log',root+'.tex'))
            for figure in figures: os.system('rm %s' % figure)
        return root+'.ps'

    def Compile2PDF(self,LatexFile,cleanafter=False,figures=[]):
        """Compiles a Latex file"""
        # IMPORT STUFF
        import os
        import string
        # END IMPORT
        
        root = LatexFile[0:string.rfind(LatexFile,'.tex')]
        
        execline1 = 'latex %s' % LatexFile
        os.system(execline1)
        execline2 = 'dvipdf %s %s' % (root+'.dvi',root+'.pdf')
        os.system(execline2)
        
        if cleanafter :
            os.system('rm %s %s %s %s' % \
            (root+'.dvi',root+'.aux',root+'.log',root+'.tex'))
            for figure in figures: os.system('rm %s' % figure)
        return root+'.pdf'


    def Ps2Pdf(self,PsFile,PdfFile,cleanafter=False):
        """Converts ps to pdf using ps2pdf"""
        # IMPORT STUFF
        import os
        # END IMPORT
        os.system('ps2pdf %s %s' % (PsFile,PdfFile))
        if cleanafter : os.system('rm %s' % PsFile)
        return None
    
    
    def GhostMergePdfs(pdfList,pdfFile):
        """ """
        tempcommand = 'ghostscript -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=%s %s'
        pdfListLine = st.join(pdfList,' ')
        os.system(tempcommand % (pdfFile,pdfListLine))
    
    def MergePdfs(self,pdfList,pdfFile):
        """This script merges a list of .pdf files to a single .pdf file. 
	It uses pdftk."""
        # IMPORT STUFF
        import os
        import sys
        #import copy
        from time import time
        # END IMPORT

        if len(pdfList) < 2:
            sys.exit('At least 2 files should be merged.\n What are you thinking on?')
    
    
        i = 1
        
        tmptemplate = 'tmp%s_latex_MergePdfs_%f.pdf' % ('%s',time())
    
        tmp2name = tmptemplate % i
        os.system('pdftk %s %s cat output %s' % (pdfList[0],pdfList[1],\
            tmp2name))
    
    
        if len(pdfList) >= 3:
    
            for i in range(2,len(pdfList)):
                tmp1name = tmptemplate% (i-1)
                tmp2name = tmptemplate % i
            os.system('pdftk %s %s cat output %s' % (tmp1name,pdfList[i],\
            tmp2name))
            os.system('rm %s' % tmp1name)
    
        last_merged = tmp2name
     
        os.system('mv %s %s' % (last_merged,pdfFile))
        os.system('rm %s' % (tmptemplate % '*',))
        
        return None
    
