.. _installation:

Installation
============


Cloning *vison* from the repository using *git*
-----------------------------------------------

If you don't have *git* installed in your system, please follow this 
`link <https://www.atlassian.com/git/tutorials/install-git>`_ first.


Here we will follow these `instructions <https://help.github.com/articles/cloning-a-repository/ for a linux system>`_ 
to clone the repository to your own computer. Follow the link for instructions in other operative systems.

Step-by-step:

* Go to https://github.com/ruymanengithub/vison.
* Click on the green "Clone or download" button.
* In the Clone with HTTPs section, click  to copy the clone URL for the repository.
* Open a Terminal.
* Change the current working directory to the location where you want the cloned directory 
  to be made.
* Type ``git clone``, and then paste the URL you copied in Step 1.
  ::
    ~$ git clone https://github.com/ruymanengithub/vison
    
* Press Enter. Your local clone will be created.

.. figure:: figs/repo_screenshot.png
    :align: center


Installation
------------

We recommend installing the code through a `conda` environment, with a specific list
of packages, so you can be sure you have all the needed dependencies.

First, if you don't have `conda` already installed in your system already, 
follow the instructions in this `link <https://conda.io/docs/install/quick.html>`_.


Installing conda and creating `vison` environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once you have succesfully installed conda, we will create an environment that will
allow you to install the pipeline and meet all its dependencies (save for SAO DS9, which 
is only used in real-time monitoring, optionally).


Step-by-Step:

* change directory to your copy of the vison repository:
  ::
    ~$ cd vison

* Under the 'conda' sub-folder, you will find several text files: 
  ::
    ~$ cd conda
    ~$ ls
    env-conda_vison_linux.txt  env-pip_vison.txt

* Then execute the following command to create a new conda environment, `vison`.
Use the OS version that may correspond in your case (by now, only linux-64 bits version available).

  ::
    ~$ conda create -n vison --file env-conda_[OS].txt

* When prompted, type "y" and return to install the listed packages.
* Activate the new environment
  ::
    ~$ source activate vison

* Update pipe
  ::
    ~$ pip install --upgrade pip

* Install the packages that are accessed via `pip`, within the conda environment:
  ::
    ~$ pip install -r env-pip_vison.txt
        

Installing `vison`
^^^^^^^^^^^^^^^^^^

Finally, to install the `vison` pipeline itself, we will go back to the folder we downloaded
from the github repository:
::
    ~$ cd ../
    ~$ ls
    conda  docs  LICENSE  manual_vison.pdf  README.md  setup.cfg  setup_distutils.py  setup.py  vison

Then do the actual installation, via:

~$ python setup.py install

Now the vison package will be accessible from anywhere in your system, whenever you start python from 
within the `vison` conda environment. For example:

* open a new terminal and go to your home directory
  ::
    ~$ cd 

* activate the vison environment:
  ::
    ~$ source activate vison

* start the python interpreter and import vison:
  ::
    ~$ source activate vison
    ~$ python
    >>> import vison
    >>> dir(vison)
    ['Eyegore', 'FlatFielding', 'Pipe', 'Report', '__all__', '__builtins__', '__doc__', '__file__', 
    '__name__', '__package__', '__path__', '__version__', '_version', 'analysis', 'blocks', 'dark', 
    'data', 'datamodel', 'eyegore', 'flat', 'image', 'inject', 'matplotlib', 'ogse', 'ogse_profiles', 
    'other', 'pipe', 'plot', 'point', 'pump', 'stop', 'support']

