Parameterization toolkits
=========================

Assigning atomic partial charges using some flavor of `AM1-BCC <https://docs.eyesopen.com/toolkits/python/quacpactk/molchargetheory.html#am1bcc-charges>`_
with `polymerist` also requires installation of some supplementary toolkits.

One can mix-and-match installing any combination of the toolkits below to taste. 
If impatient or indifferent, one might opt for a "shotgun" approach and install all 3 at once with the following commands:

.. code-block:: console

    mamba activate polymerist-env
    mamba install -c openeye openeye-toolkits
    mamba install -c conda-forge espaloma_charge "torchdata<=0.9.0"
    mamba install -c conda-forge openff-nagl "torchdata<=0.9.0"

`OpenEye toolkits <https://docs.eyesopen.com/toolkits/python/intro.html>`_ (closed-source)
------------------------------------------------------------------------------------------
These toolkits are required to perform explicit 
`AM1-BCC charge assignment with conformer selection <https://docs.eyesopen.com/toolkits/python/quacpactk/molchargetheory.html#elf-conformer-selection>`_
and enhance OpenFF's conformer generation, but **`polymerist` does not *require* OpenEye to work**. 

If you already have (or can obtain) an `OpenEye license <https://docs.eyesopen.com/toolkits/python/quickstart-python/license.html>`__
and wish to install the OpenEye toolkits individually, follow the 
`install instructions provided by OpenEye <https://docs.eyesopen.com/toolkits/python/quickstart-python/linuxosx_anaconda.html#:~:text=Install%20the%20OpenEye%20Python%20Toolkits%20into%20the%20new%20environment%3A>`_.

`Espaloma-charge <https://github.com/choderalab/espaloma-charge>`_
------------------------------------------------------------------
This is a standalone graph neural network (GNN) model [`Wang et. al. <https://pubs.acs.org/doi/10.1021/acs.jpca.4c01287>`_]
which can assign atomic partial charges trained on AM1-BCC data extremely rapidly. To install individually, follow the
`install instructions provided by Espaloma's developers <https://github.com/choderalab/espaloma-charge?tab=readme-ov-file#installation>`_.

`OpenFF NAGL <https://docs.openforcefield.org/projects/nagl/en/latest/index.html>`_
-----------------------------------------------------------------------------------
This is an OpenFF-specific GNN based on similar architecture to Espaloma with a generally
better validated partial charge model. To install individually, follow the
`install instructions provided by the OpenFF NAGL documentation <https://docs.openforcefield.org/projects/nagl/en/latest/installation.html#:~:text=If%20you%20prefer%2C%20NAGL%20may%20be%20installed%20into%20the%20current%20environment%3A>`_.
