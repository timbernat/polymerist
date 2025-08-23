Installing from source
======================

Polymerist can also be installed directly from its source code in the
`polymerist GitHub repository <https://github.com/timbernat/polymerist>`_.
To install, execute the following set of terminal commands in whichever
directory you'd like the installation to live on your local machine:

Mamba installation from source
------------------------------

.. code-block:: console

    git clone https://github.com/timbernat/polymerist
    cd polymerist
    mamba env create -n polymerist-env -f devtools/conda-envs/release-build.yml
    mamba activate polymerist-env
    pip install .

Once the source install is complete, you no longer need the clone
of the polymerist repo and can remove it from your file system.
