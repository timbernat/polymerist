Polymer-Oriented LibrarY of Monomer Expression Rules and In-silico Synthesis Tools
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/timbernat/polymerist/workflows/CI/badge.svg)](https://github.com/timbernat/polymerist/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/timbernat/polymerist/main/graph/badge.svg)](https://codecov.io/gh/timbernat/polymerist/branch/main)

A unified set of tools for monomer template generation, topology building, force-field parameterization, and MD simulations of general organic polymer systems within the OpenFF framework

Source code for [Davel, Connor M., Bernat, Timotej, Wagner, Jeffrey R., and Shirts, Michael R., "Parameterization of General Organic Polymers within the Open Force Field Framework"](https://pubs.acs.org/doi/10.1021/acs.jcim.3c01691)

![abstract](docs/_static/polymer_param_graphic_TOC.png)


## Installation

Currently, this package only supports a "dirty" developer install using conda/mamba and `pip`.  
First, you will need to install the `conda` environment manager (either the lightweight [Miniconda Distribution](https://docs.anaconda.com/free/miniconda/miniconda-install/) (recommended) or the bulkier [Anaconda Distribution](https://www.anaconda.com/download)) if you don't already have it.
Further, it is highly recommended (but optional) that you also have the `mamba` package manager installed (either through [Miniforge or Conda](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html)); this will greatly accelerate download times.  

Once `mamba` is installed, you can proceed with the dirty `polymerist` install into a safe virtual environment (named "polymerist-env" here).
To install, execute the following set of commands in a command line interface (CLI) in whichever directory you'd like the dev installation to live:
```sh
git clone https://github.com/timbernat/polymerist
cd polymerist
mamba env create -n polymerist-env -f devtools/conda-envs/polymerist-env.yml
mamba activate polymerist-env
pip install -e .
```
The third command will take **_at least_** a few minutes, and will make the CLI terminal quite busy; remain calm, that's normal!

Equivalent commands using just `conda` (in case `mamba` has not been installed) are below. These will perform the same installation, just much more slowly:
```sh
git clone https://github.com/timbernat/polymerist
cd polymerist
conda env create -n polymerist-env -f devtools/conda-envs/polymerist-env.yml
conda activate polymerist-env
pip install -e .
```
As an optional last step, it is recommended that you correctly set up paths to your [OpenEye License](https://docs.eyesopen.com/toolkits/python/quickstart-python/license.html), if you have access to one.
Portions of conformer-generation and partial-charge assignment in `polymerist` will work more effectively with the OpenEye toolkit installed and licensed, but 'polymerist' is set up to not require these closed-source dependencies.

From here, you should be able to run `polymerist`-dependent scripts in the polymerist-env virtual environment active, either from the command line or from a [Jupyter Notebook](https://jupyter.org/)!


### Copyright

Copyright (c) 2024, Timotej Bernat (timotej.bernat@colorado.edu)


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
