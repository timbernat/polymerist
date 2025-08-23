Prerequisites
=============

Python
------
Currently, polymerist is only compatible `Python 3.11 <https://www.python.org/downloads/release/python-3110/>`_;
**attemping to install with any other version of Python will not work.**

OS
--
polymerist is compatible with Linux (recommended) and Mac machines capable of installing Python 3.11.
Due to `lack of support from AmberTools <https://ambermd.org/InstWindows.php>`_, direct installation on Windows machines is not supported;
However, this can easily be circumvented by using the `Windows Subsystem for Linux (WSL2) <https://learn.microsoft.com/en-us/windows/wsl/install>`_

Python package manager
----------------------
Before proceeding with installation, ensure you have some iteration of a
Python package and environment management system installed on your machine.
If you don't already have one installed, it's recommended you install either of:

* mamba through `Miniforge <https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html>`_
* conda either as the lightweight `Miniconda <https://docs.anaconda.com/free/miniconda/miniconda-install/>`_ or the significantly bulkier `Anaconda <https://www.anaconda.com/download>`_.

mamba is very rapid and is the recommended package manager for this install; **if you opt to
use `conda` over `mamba`, be prepared for a markedly slower and more tedious install process!**
Users with a pre-existing conda installation can `still install mamba <https://anaconda.org/conda-forge/mamba>`_, so there's really no reason not to use it.

Once you have a package manager installed, you may proceed with one of the installation methods detailed below.