Getting Started
===============

This page details how to get started with Polymer-Oriented LibrarY of Monomer Expression Rules and In-silico Synthesis Tools. 

Installation
############
Polymerist has the following dependencies:

* Python 3.11.z (no later version supported)
* NumPy
* SciPy
* RDKit

To install, run the following commands:
::

    git clone https://github.com/timbernat/polymerist
    cd polymerist
    pip install .

Usage
#####
.. code-block:: python
    
    import polymerist as ps

    print(ps.pascal(5))