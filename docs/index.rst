.. vison documentation master file, created by
   sphinx-quickstart on Thu Feb  2 10:44:26 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to vison's documentation!
=================================

Contents:

.. toctree::
   :maxdepth: 2
   
   source/includeme

:Author: Ruyman Azzollini
:Contact: r.azzollini_at_ucl.ac.uk
:issue: |release|
:version: |version|
:date: |today|

This Python package "vison" is the software pipeline that has been used at MSSL for ground calibration of the VIS detection chains (12 + 2 spares), including one ROE, one RPSU and three CCDs each.

Installation
============

The package is distributed via github. The repository is hosted at:

https://github.com/ruymanengithub/vison


Detailed instructions:

.. toctree::
    :maxdepth: 4
    
    installation

Dependencies
------------

Instructions to acquire a copy of the "conda" environment that provides almost all 
dependencies is included in the package. See :ref:`installation` instructions for 
details. The only package that will not be installed by this means is SAO-DS9,
which must be installed separetly, in order to be able to use some of the
interactive inspection capabilities of Eyegore (data acquisition monitoring).

How to Use it
=============

Some kind of **cook-book**. 

*A bit like when your grandmother tells you to put a pinch of sugar and heat it for a while... that kind of accuracy.*

.. toctree::
    :maxdepth: 4

    cookbook

Guide through the code
======================

In this section we provide an overal description of the capabilities and organisation of the code.

.. toctree::
    :maxdepth: 4
    
    guide


Pipeline Core
=============

Pipeline master classes.

.. toctree::
    :maxdepth: 4
    
    pipe


Data Model
==========

Modules with classes to hold data model for inputs and outputs: exposure log,
HK files, FITS files, etc.

.. toctree::
   :maxdepth: 4
   
   datamodel


Analysis (Shared)
=================

.. toctree::
    :maxdepth: 4

    analysis


Charge Injection Tools
======================

.. toctree::
    :maxdepth: 4
    
    inject


"Flat" Acq. Analysis Tools
==========================

.. toctree::
    :maxdepth: 4
    
    flat


FPA Tests
=========

.. toctree::
    :maxdepth: 4
    
    fpatests

Image
=====

.. toctree::
    :maxdepth: 4

    image

Metatests
=========

.. toctree::
    :maxdepth: 4
    
    metatests

Monitoring ("Eyegore")
======================

Tools to monitor data acquisition on real time: plots of HK, auto-updating of
visual display of Exposure Log with some interactive capabilities, and display of latest images.

.. toctree::
    :maxdepth: 4

    eyegore

OGSE
====

OGSE stands for Optical Ground Support Equipment.

.. toctree::
    :maxdepth: 4

    ogse

Plotting
========

General use plotting facilities.

.. toctree::
    :maxdepth: 4

    plot

Point-Source Analysis
=====================

.. toctree::
    :maxdepth: 4

    point

Scripts
=======

These are pipeline scripts, not the Test Scripts (for those keep scrolling down).

.. toctree::
   :maxdepth: 4
   
   scripts


Support Code
============

.. toctree::
   :maxdepth: 4
   
   support

Unit Testing
============

.. toctree::
    :maxdepth: 4

    utests



Test Scripts
============

These are the scripts that hold the description, execution, data validation and
analysis of the tests that make the campaign. They are served by the
infrasctructure and tools provided by the pipeline.


Charge Injection Scripts
------------------------

.. toctree::
    :maxdepth: 4

    chinj_scripts


Dark Scripts
--------------

.. toctree::
    :maxdepth: 4
    
    dark_scripts


Flat-Illumination Scripts
-------------------------


.. toctree::
    :maxdepth: 4

    flat_scripts


Point-Source Scripts
--------------------

.. toctree::
    :maxdepth: 4

    point_scripts

Trap-Pumping Scripts
--------------------

.. toctree::
    :maxdepth: 4
    
    tp_scripts



Other Test Scripts
------------------

.. toctree::

    other_scripts






Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


