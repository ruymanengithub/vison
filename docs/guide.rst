.. _guide:


Code Guide
==========

Here we provide an overal description of the code, what it does and how it is organised, to help the user / developer understand how things work.

Capabilities of the Code
------------------------

First, it is convenient to know what the code here provided can do:

1. Write **excel scripts for the acquisition of data** using ELVIS. Each test has an associated class, and this class has a method which generates the description of the test that is then converted to excel. The test description can be modified at runtime by means of test object parameters. These test descriptions are a key element of the pipeline. Monitoring and analysis capabilities depend on them (i.e. knowing what *should* happen).
2. **Monitor the data acquisition** of the calibration campaign in real time. This is done through the sub-package "eyegore". It can monitor the data flow, values in the HK data stream (cheking value against limits), display images and the exposure log, and issue warnings via sms and e-mail.
3. **Analyse the data** from the tests. This is the core functionality of the pipeline.
4. Produce calibration / analysis results in a range of formats: **Json, FITS, excel**, etc.
5. Produce **test reports** for each test in pdf (via LaTeX).
6. Collate and produce summary reports from a set of tests across all calibrated blocks (**metatests**).
7. Analyse data from the whole VIS FPA (**fpatests**).
8. **Simulate data** (very basic capabilities).


Code Architecture
-----------------

**TBW**


Code Flow
---------

It is perhaps easier to describe what the pipeline does, and how it is organised, following what it does when we try to analyse a data-set.




Data Models
-----------

See :ref:`datamodel`


Reusing / adapting the code
---------------------------