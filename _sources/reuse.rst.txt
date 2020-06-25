
.. _reuse:


Repurposing the vison package
=============================


Could we reuse this package to write scripts, monitor the data acquisition, analyze data, etc. from other instruments? Probably, but it is not straight-forward, at best. Let me first describe the current situation, enumerating all aspects that make the code specific to the VGCC:

* The test definitions are specific to the VGCC, obviously. These specificity is mostly embodied in the methods of the Task classes (specifically through the .build_scriptdict() method, but also in the other methods which expect a particular test structure). Then, the campaign description itself, with the specific values of free parameters of the tests (e.g. exposure times, or number of frames, as in some tests), is also specific. This second aspect is much easier to modify than the test/task classes, of course.
* The main data structure for each task is the *dd* attribute (DataDict instance). This comes shaped by the EXPLOG output by ELVIS. The specific column names and contents in the EXPLOG, and the HK files can be relatively easy modified, as their format and contents are handled through specific modules/functions, but probably quite a few things could break if the format of the EXPLOG or HK files was different (if they are not there at all, then that probably would make the possibility of "migration" to vanish... it could be better to start from scratch).
* Insisting on the DataDict catalog / data structure, this is a *monster* of its own, and I was to write the code again, I'd consider using something less contrived. Probably I'd go for a pandas or an astropy table structure, or an object that would encapsulate an object of that kind. It could be done, and it could be worth it, if the code was to be reused for more than a single mission, but that would require knowledge of the code, not to waste too much time finding your way through the code.
* The analysis of the images relies mostly on an object, the CCD class, which is shaped according to the Euclid CCDs. It would not be too difficult to repurpose this as a CCD with other dimensions and architecture, because things are mostly defined just once, but it's not trivial either.
* The fact that there are 3 CCDs per ROE in the VIS detection chain is probably hardwired in more places than it should be, and may pose a problem when migrating to a system with a different number of CCDs per ROE / detection chain.
* Plotting is done through figure objects, and these are shaped for the specific case of EUCLID. As an example, one of the most repeated plots is shaped to have 12 panels, as there are 3 CCDs per beam, 4 quadrants per CCD.


Other considerations:

* The code is still poorly documented.
* The code works for what it was intended, but it wasn't designed as a general purpose code. This means that in some cases, decisions were made taking into account limited resources (time), and in many cases, the working solution was preferred over the general solution.


You are still not deterred and want to reuse the code
-----------------------------------------------------

Brave you! Well, here is some kind of to-do list to accomplish that:

* Collect all the information regarding the description of the detection chains, the OGSE, and the outputs from the lab acquisition software. Maket it easy to share and edit collaboratively. A google drive could help you with this, just sayin'...
* Clone the vison project to another project in github, or wherever you want to host the new project (version control is a must, and having a remote repository is an added safety measure against data loss).
* **No, it's not good idea to split vison in a EUCLID area and some additional sub-packages for other(s) mission(s)**. It will get far too messy in no time. Let *vison* die, let it have (improved, more versatile) offspring. It's absolutely recommended to start a new project.
* First, I'd tackle the interfaces to the new? lab output files, starting with the EXPLOG and the HK. Adapt the existing modules to ingest those new files, test that it works.
* Then go for the CCDs, and create a class "XCCD" that inherits from CCD, and try to make it work for the new CCD images just by adjusting the values of pre-scan, over-scan, quadrant lists, etc. If that doesn't suffice, keep working until you get a CCD class that works for you, and you can at least do the same analysis that was possible using CCD for Euclid (get statitistics from selected regions, get 2D models, do some arithmetics with the images, collapse them along different axes, etc.).
* Time to create the first Task class. Consider the easiest test, ideally one that may be common with Euclid. For example a bias acquisition test? Try to create a class for your specific test, and make it write the acquisition script for your test.
    * Modify the build_scriptdict() method, and all the other things you may need to accomplish this.
* Then, you could create a second class, for a different test (flat-fields?).
* Now you could try to create a minimum campaign with just these 2 tests.
* If you have available already acquired test data, try to run the pipeline to ingest that campaign and use the newly created test classes, and do the analysis... good luck!
* If you don't yet have the data, you could try to generate it yourself, with the data simulator (which may also need adaptation, though after modifying the other modules, this will be easy cake).


I just want to write acquisition scripts for "ELVIS"
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Then use Pandas, it can write excel files, for example. 

If you want to also have classes that hold the description of the test, produce the script according to this description, and then you want to analyse the data, then, well, you'll have to follow the path above, or something similar (if reusing vison, that is, of course).


