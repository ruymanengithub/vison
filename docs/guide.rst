.. _guide:


Code Guide
==========

Here we provide a succint description of the code, what it does and how it is organised, to help the user / developer understand how things work, to better make use of it / repurpose it.

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

Here it is convenient to introduce some nomenclature regarding the executiong of the pipeline that you may see used in the naming of classes and functions therein:

* A **Task** is basically a test (e.g. BIAS02), with a description of how the test data should be acquired, methods to analyse the data once acquired, and others to plot, produce reports, check compliances, etc. 
* The execution of a Task is broken down in subtasks, which are methods of the the Task classes.
* A pipeline is a sequential execution of a number of a tasks. In the pipeline it takes shape as a class (Pipe) that has to be instantiated with inputs to perform analysis on a number of tests, using the data stored somewhere.


Code Flow
---------

It is perhaps easier to describe what the pipeline does, and how it is organised, following what it does when we use to perform different tasks.

Data-set Analysis
^^^^^^^^^^^^^^^^^

Let's first try to analyse a data-set.

We will call the script **vison_run** which instantiates a Pipeline object, loads it with the tests that we are going to process, the inputs to the tasks for those tests, and then runs the pipeline object. Let's go step by step with an example.

::
    ~$ vison_run -y [vison_config.py] -R [SESSION] -l -t [TAG]

Here vison_config.py stands for a python script with inputs (more on that soon), SESSION is a name to select the acquisition session within the the configuration file we want to select for analysis (there usually are several sessions within a configuration script, as there are data acquisition sessions in a multi-day campaign), and TAG is just a character string to label the directory with results, and the text log file, for ease of identification.

Before we go on, some basic notions regarding the organisation of the GCC:

* the campaign was sub-divided in campaigns for each block.
* within each block-campaign there were sessions/runs (-R comes from "run") in which several tests are executed one after another, autonomously by ELVIS.

For each block-campaign we created a vison_config script, with the name of the block in question and the date of the campaign (e.g. vison_config_EINSTEIN_JUL19.py). This script is important because it is required to run the pipeline, and because it serves as registry of the inputs used to run it.

Let's have a look at the contents of a vison_config.py script, in simplified form.

.. literalinclude:: vison_config_example.py
   :language: python

When vison_run is executed, this literally *executes* the vison_config script which generates a dictionary called **inputdict**. This dictionary has all the information necessary for the pipeline to run:

* in which facility (*chamber*) was the data acquired (relevant to know exposure times, for example).
* what tasks (tests) are to be executed/analysed.
* what values to assign to the free parameters of the tasks.
* what hardware are we testing.
* where to find the input data, and where to put the output results.
* where to find input calibration data products (e.g. cosmetics masks).

When vison_config.py is executed, it starts from:

::
    if __name__ == '__main__'

From there, it calls the function *get_config_dict* to create a standard version of the inputs dictionary, setting up the tasks and their standard inputs. Then it calls *add_RUN_specifics* to add data locations, OBSID ranges, and apply a selector of sub-tasks to execute for each task.

Going back to *vison_run*, once the inputs in the configuration file are ingested, the next important thing is to create an instance of the pipeline class:

::
    pipe = Pipe(inputdict, dolog=dolog, drill=drill,
                    debug=debug, startobsid=startobsid,
                    processes=multithread, tag=tag,
                    cleanafter=cleanafter)

This takes as input the input dictionary (inputdict) and other keywords to control the execution of the pipeline.

Then, the pipeline is executed, either in wait-for-data mode (in parallel with acquisition), or directly (assuming all data has already been acquired previously):

::
    if wait:
        pipe.wait_and_run(dayfolder, elvis=elvis)
    else:
        pipe.run()




Data Models
-----------

See :ref:`datamodel`


Reusing / adapting the code
---------------------------