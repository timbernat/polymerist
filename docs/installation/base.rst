Base installation
=================

A fully-featured install of polymerist into a fresh virtual environment
(called "polymerist-env" here, but can be renamed to whatever you like)
can be obtained by running the following shell commands:

.. code-block:: bash 

    mamba create -n polymerist-env python=3.11
    mamba activate polymerist-env
    pip install polymerist
    mamba install -c conda-forge openff-toolkit mbuild openbabel "packmol<=20.15.1"

The final mamba install is done as a one-liner (rather than in parts)
to reconcile any discrepancies between packmol versions between the
3 conda-only packages required (namely the openff-toolkit, mbuild, and openbabel).

The packmol pin is to ensure solvation utilities work correctly; this will save you headaches,
any later versions give inscrutable PDB errors where there previously were none.

Installing the 'openff-toolkit <https://github.com/openforcefield/openff-toolkit`_ will take
**_at least_** a few minutes, and will make the terminal output quite busy; remain calm, that's normal! 

Test installation
-----------------
To test that the installation worked, you can trying to run the following snippet in python:

.. code-block:: console
    
    mamba activate polymerist-env; python
    >>> import polymerist as ps
    >>> print(ps.pascal(5))
        1    
       1 1   
      1 2 1  
     1 3 3 1 
    1 4 6 4 1

Extended install
----------------
An extended install with `Jupyter Notebook <https://jupyter.org/>`_ support, molecular visualization
capability, and chemical database querying capability can be obtained with similar commands.
Testing of installation is analogous to the python code above.

.. code-block:: bash 

    mamba create -n polymerist-env python=3.11
    mamba activate polymerist-env
    pip install polymerist[interactive,chemdb]
    mamba install -c conda-forge openff-toolkit mbuild openbabel "packmol<=20.15.1"

Conda install (**NOT** recommended)
-----------------------------------
As mentioned, the openff-toolkit install above is pretty slow with mamba; with conda, it **will** be glacial
Nevertheless, we provide equivalent conda installation instruction here for users who are too stubborn to use mamba instead.
These will perform the same installation, just much more slowly:
:: 

    conda create -n polymerist-env python=3.11
    conda activate polymerist-env
    pip install polymerist[interactive,chemdb]
    conda install -c conda-forge openff-toolkit mbuild openbabel "packmol<=20.15.1"