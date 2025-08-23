========
Examples
========

Tutorial examples are provided as `Jupyter <https://jupyter.org/>`__ notebooks in a separate
`polymerist-examples <https://github.com/timbernat/polymerist_examples>`__ repository.
These contain a series of ready-to-execute notebook accompanied by explanations of the
"how" and "why" of running workflows with the core parts of polymerist.

The goal of the examples is to show you how to set up and run molecular dynamics
simulations of polymer systems starting from nothing but knowledge of monomers,
polymerization mechanism, and any co-molecules (e.g. solvent) you wish to model.

1 - `Building polymers <https://github.com/timbernat/polymerist_examples/blob/main/1-polymerization/1.0-index.ipynb>`_
======================================================================================================================
* 1.1 - `Introduction: nylon polymers <https://github.com/timbernat/polymerist_examples/blob/main/1-polymerization/1.1-nylon_basics.ipynb>`_
* 1.2 - `Simple vinyl polymers and autopolymerization <https://github.com/timbernat/polymerist_examples/blob/main/1-polymerization/1.2-vinyl_autopolymerization.ipynb>`_
* 1.3 - `Kapton polyimides and with multiple intermonomer bonds <https://github.com/timbernat/polymerist_examples/blob/main/1-polymerization/1.3-polyimide_multibond_cycles.ipynb>`_
* 1.4 - `Crosslinkable MPD-TMC polyamides <https://github.com/timbernat/polymerist_examples/blob/main/1-polymerization/1.4-MPD-TMC_polyamides.ipynb>`_
* 1.5 - `PEG-PLGA block copolymers <https://github.com/timbernat/polymerist_examples/blob/main/1-polymerization/1.5-PEG-PLGA_copolymers.ipynb>`_
* 1.6 - `conjugated thiophenyl polymers with arbitrary sidechains <https://github.com/timbernat/polymerist_examples/blob/main/1-polymerization/1.6-functionalized_polythiophenes.ipynb>`_

2 - `Preparing systems containing polymers <https://github.com/timbernat/polymerist_examples/blob/main/2-preparation/2.0-index.ipynb>`_
=======================================================================================================================================
* 2.1 - `Loading a polymer structure into OpenFF <https://github.com/timbernat/polymerist_examples/blob/main/2-preparation/2.1-loading_polymer_topology.ipynb>`_
* 2.2 - `Polymer metadata and partial charge assignment <https://github.com/timbernat/polymerist_examples/blob/main/2-preparation/2.2-preparing_individual_polymers.ipynb>`_
* 2.3 - `Solvation and packing of polymer melts <https://github.com/timbernat/polymerist_examples/blob/main/2-preparation/2.3-melt_packing_and_solvation.ipynb>`_
* 2.4 - `Reduction Charge Transfer (RCT) for generating custom library charges <https://github.com/timbernat/polymerist_examples/blob/main/2-preparation/2.4-RCT_demo.ipynb>`_

3 - `Running polymer simulations <https://github.com/timbernat/polymerist_examples/blob/main/3-simulation/3.0-index.ipynb>`_
============================================================================================================================
* 3.1 - `Exporting polymer systems to common MD engines <https://github.com/timbernat/polymerist_examples/blob/main/3-simulation/3.1-MD_export_with_Interchange.ipynb>`_
* 3.2 - `Reproducibly serializing OpenMM simulations <https://github.com/timbernat/polymerist_examples/blob/main/3-simulation/3.2-serializable_simulation_parameters.ipynb>`_
* 3.3 - `Running series of OpenMM simulations <https://github.com/timbernat/polymerist_examples/blob/main/3-simulation/3.3-running_openmm_simulations.ipynb>`_
* 3.4 - `A start-to-finish polymer simulation workflow <https://github.com/timbernat/polymerist_examples/blob/main/3-simulation/3.4-full_workflow_demo.ipynb>`_

`Miscellaneous <https://github.com/timbernat/polymerist_examples/blob/main/4-misc/4.0-index.ipynb>`_
====================================================================================================
* `Exporting arbitrary Python objects to JSON <https://github.com/timbernat/polymerist_examples/blob/main/4-miscellanous/jsonification.ipynb>`_
* `Automated ring piercing detection (PINPRICS) <https://github.com/timbernat/polymerist_examples/blob/main/4-miscellanous/ring-piercing.ipynb>`_
* `Loading boron-containing molecules into OpenFF <https://github.com/timbernat/polymerist_examples/blob/main/4-miscellanous/openff_boron_demo.py>`_