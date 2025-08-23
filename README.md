Polymer-Oriented LibrarY of Monomer Expression Rules and In-silico Synthesis Tools
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/timbernat/polymerist/workflows/CI/badge.svg)](https://github.com/timbernat/polymerist/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/timbernat/polymerist/main/graph/badge.svg)](https://codecov.io/gh/timbernat/polymerist/branch/main)

A unified set of tools for setting up molecular dynamics simulations of general organic polymer systems. Based upon concepts introduced in ["Parameterization of General Organic Polymers within the Open Force Field Framework" (Davel, Connor M., Bernat, Timotej, Wagner, Jeffrey R., and Shirts, Michael R.)](https://pubs.acs.org/doi/10.1021/acs.jcim.3c01691)

![abstract](docs/_static/polymer_param_graphic_TOC.png)

## Features
Includes functionality for:
* Generating chemical information-rich monomer residue templates
* Enumeration of all possible repeat units of a polymer ensemble given initial monomers and a target polymerization mechanism
* Building of linear homopolymers and copolymers (both topologies and coordinates)
* Solvating polymer systems
* Packing of many polymer chains into melt-like boxes
* Force-field parameterization within the [OpenFF](https://openforcefield.org/about/organization/) framework
* Interfaces to semi-empirical and graph neural network-based atomic partial charge assignment
* Reproducible simulation specification API
* Interfaces to OpenMM for running serial simulations with different thermodynamic parameters
* Much more!

## Documentation
Complete documentation for polymerist can be found on the [`polymerist` ReadTheDocs page](https://polymerist.readthedocs.io/en/docs/)
  
## Examples
Examples of how to import and invoke the core features of `polymerist` can be found in the accompanying [polymerist_examples repository](https://github.com/timbernat/polymerist_examples).

## Requirements
### OS
`polymerist` is compatible with Linux (recommended) and Mac machines capable of installing Python 3.11. Due to [lack of support from AmberTools](https://ambermd.org/InstWindows.php), direct installation on Windows machines is not supported; **however**, this can easily be circumvented by using the [Windows Subsystem for Linux (WSL2)](https://learn.microsoft.com/en-us/windows/wsl/install)

### Python package manager
Before proceeding with installation, ensure you have some iteration of a Python package and environment management system installed on your machine. We recommend using `mamba`, installable via [Miniforge](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html). [Miniconda](https://docs.anaconda.com/free/miniconda/miniconda-install/) and other [Anaconda](https://www.anaconda.com/download) manager also work; however **if you opt to use `conda` over `mamba`, be prepared for a markedly slower and more tedious install process!**

## Installation
The distribution for `polymerist` is hosted on [PyPI](https://pypi.org/project/polymerist/). For detailed instruction on how to install polymerist and other required toolkits, see [the installation docs](https://polymerist.readthedocs.io/en/docs/installation/index.html). Here we provide an abridged summary to get you going

### Base install
A fully-featured install in a safe virtual environment (named "polymerist-env", here) can be obtained by running the following terminal commands:
#### Mamba install (basic)
```bash
mamba create -n polymerist-env python=3.11
mamba activate polymerist-env
pip install polymerist
mamba install -c conda-forge openff-toolkit mbuild openbabel "packmol<=20.15.1"
```

#### Mamba install (extended)
An extended install with [Jupyter Notebook](https://jupyter.org/) support, molecular visualization capability, and chemical data querying capability can be obtained very similarly:
```bash
mamba create -n polymerist-env python=3.11
mamba activate polymerist-env
pip install polymerist[interactive,chemdb]
mamba install -c conda-forge openff-toolkit mbuild openbabel "packmol<=20.15.1"
```

### Testing installation
To see if the installation was successful, one can run the following short set of commands which should yield the outputs shown:
```python
mamba activate polymerist-env; python
>>> import polymerist as ps
>>> print(ps.pascal(5))
    1    
   1 1   
  1 2 1  
 1 3 3 1 
1 4 6 4 1
```

### Parameterization toolkits
Assigning atomic partial charges using some flavor of [AM1-BCC](https://docs.eyesopen.com/toolkits/python/quacpactk/molchargetheory.html#am1bcc-charges) and enhanced conformer generation also requires installation of some supplementary toolkits. These can be installed as follows:
```bash
mamba activate polymerist-env
mamba install -c openeye openeye-toolkits
mamba install -c conda-forge espaloma_charge "torchdata<=0.9.0"
mamba install -c conda-forge openff-nagl "torchdata<=0.9.0"
```

### Installing from source (optional)
`polymerist` and all required dependencies can also be installed directly from the source code in this repository.
To install, execute the following set of terminal commands in whichever directory you'd like the installation to live on your local machine:
```bash
git clone https://github.com/timbernat/polymerist
cd polymerist
mamba env create -n polymerist-env -f devtools/conda-envs/release-build.yml
mamba activate polymerist-env
pip install .
```

### Developer installation (for advanced users only)
Those developing for `polymerist` may like to have an editable local installation, in which they can make changes to the source code and test behavior changes in real-time. This type of installation proceeds as follows:
```bash
git clone https://github.com/timbernat/polymerist
cd polymerist
mamba env create -n polymerist-dev -f devtools/conda-envs/dev-build.yml
mamba activate polymerist-dev
pip install -e . --config-settings editable_mode=strict
```
The `--config-settings editable_mode` flag in the final line allows this installation to "play nicely" with PyLance, making auto-completion and navigation to source code much easier for VSCode users. **It is optional, and can be removed if this compatibility is not desired**


## Copyright
Copyright (c) 2024, Timotej Bernat (timotej.bernat@colorado.edu)

## Acknowledgements
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
