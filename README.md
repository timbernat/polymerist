Polymer-Oriented LibrarY of Monomer Expression Rules and In-silico Synthesis Tools
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/timbernat/polymerist/workflows/CI/badge.svg)](https://github.com/timbernat/polymerist/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/timbernat/polymerist/main/graph/badge.svg)](https://codecov.io/gh/timbernat/polymerist/branch/main)

A unified set of tools for monomer template generation, topology building, force-field parameterization, and MD simulations of general organic polymer systems within the OpenFF framework

Source code based upon concepts introduced in ["Parameterization of General Organic Polymers within the Open Force Field Framework" (Davel, Connor M., Bernat, Timotej, Wagner, Jeffrey R., and Shirts, Michael R.)](https://pubs.acs.org/doi/10.1021/acs.jcim.3c01691)

![abstract](docs/_static/polymer_param_graphic_TOC.png)


## Installation

Currently, this package only supports a "dirty" `pip` installation from a clone of this GitHub repo.
First, you will need to install some iteration of the `conda` environment manager. Possible builds (in decreasing order of recommendation) are:
- The very-rapid `mamba` package manager; installable either through [Miniforge or Conda](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html)
- The lightweight [Miniconda Distribution](https://docs.anaconda.com/free/miniconda/miniconda-install/)
- The bulky [Anaconda Distribution](https://www.anaconda.com/download)
Once you install one of these (or if you already have one installed on your system), you may proceed with on of the installation modes.

### 1) User Install

For nominal use, `polymerist` can be installed into a safe virtual environment, named "polymerist-env" here.
To install, execute the following set of commands in a command line interface (CLI) in whichever directory you'd like the dev installation to live:

#### Mamba install:
```sh
git clone https://github.com/timbernat/polymerist
cd polymerist
mamba env create -n polymerist-env -f devtools/conda-envs/release-build.yml
mamba activate polymerist-env
pip install .
```
The third command will take **_at least_** a few minutes, and will make the CLI terminal quite busy; remain calm, that's normal!

#### Vanilla Conda install:
Equivalent commands using just `conda` (in case `mamba` has not been installed) are below. These will perform the same installation, just much more slowly:
```sh
git clone https://github.com/timbernat/polymerist
cd polymerist
conda env create -n polymerist-env -f devtools/conda-envs/release-build.yml
conda activate polymerist-env
pip install .
```

### 1.1) OpenEye License
As an optional last step, it is recommended that you correctly set up paths to your [OpenEye License](https://docs.eyesopen.com/toolkits/python/quickstart-python/license.html), if you have access to one.
Portions of conformer-generation and partial-charge assignment in `polymerist` will work more effectively with the OpenEye toolkit installed and licensed, but 'polymerist' is set up to not require these closed-source dependencies.

From here, you should be able to run `polymerist`-dependent scripts in the polymerist-env virtual environment active, either from the command line or from a [Jupyter Notebook](https://jupyter.org/)!

### 2) Dev install

Those developing for `polymerist` may like to have an editable local installation, in which they can make changes to the source code and test behavior changes in real-time.
In this case, one requires an "editable build" which does **NOT** live in the site_packages dir of the created env. This type of installation proceeds as follows:
```sh
git clone https://github.com/timbernat/polymerist
cd polymerist
mamba env create -n polymerist-dev -f devtools/conda-envs/dev-build.yml
mamba activate polymerist-dev
pip install -e . --config-settings editable_mode=strict
```
The `--config-settings editable_mode` flag in the final line allows this install to "play nicely" with PyLance, making auto-completion and navigation to source code much easier for VSCode users.
It is optional, and can be removed if this compatibility is not desired.

### Copyright

Copyright (c) 2024, Timotej Bernat (timotej.bernat@colorado.edu)


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
