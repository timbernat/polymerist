Developer installation
======================

Those developing for `polymerist` may like to have an editable local installation,
in which they can make changes to the source code and test behavior changes in real-time.

In this case, one requires an "editable build" which mirrors the source files that live
in the site_packages directory of the created environment. This type of installation proceeds as follows:

.. code-block:: console
    
    git clone https://github.com/timbernat/polymerist
    cd polymerist
    mamba env create -n polymerist-dev -f devtools/conda-envs/dev-build.yml
    mamba activate polymerist-dev
    pip install -e . --config-settings editable_mode=strict

The `--config-settings editable_mode` flag in the final line allows this installation to 
"play nicely" with PyLance, making auto-completion and navigation to source code much easier
for VSCode users. **It is optional, and can be removed if this compatibility is not desired**