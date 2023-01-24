============================
SEAMM Read Structure Plug-in
============================

.. image:: https://img.shields.io/github/issues-pr-raw/molssi-seamm/read_structure_step
   :target: https://github.com/molssi-seamm/read_structure_step/pulls
   :alt: GitHub pull requests

.. image:: https://github.com/molssi-seamm/read_structure_step/workflows/CI/badge.svg
   :target: https://github.com/molssi-seamm/read_structure_step/actions
   :alt: Build Status

.. image:: https://codecov.io/gh/molssi-seamm/read_structure_step/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/molssi-seamm/read_structure_step
   :alt: Code Coverage

.. image:: https://github.com/molssi-seamm/quickmin_step/workflows/CodeQL/badge.svg
   :target: https://github.com/molssi-seamm/quickmin_step/security/code-scanning
   :alt: Code Quality

.. image:: https://github.com/molssi-seamm/read_structure_step/workflows/Release/badge.svg
   :target: https://molssi-seamm.github.io/read_structure_step/index.html
   :alt: Documentation Status

.. image:: https://pyup.io/repos/github/molssi-seamm/read_structure_step/shield.svg
   :target: https://pyup.io/repos/github/molssi-seamm/read_structure_step/
   :alt: Updates for Dependencies

.. image:: https://img.shields.io/pypi/v/read_structure_step.svg
   :target: https://pypi.python.org/pypi/read_structure_step
   :alt: PyPi VERSION

A SEAMM plug-in to read common formats in computational chemistry

The current version uses Open Babel as an engine to parse various
formats such as PDB, Mol2 or XYZ and transform them into the SEAMM
structure format for further use in SEAMM flowcharts.

* Free software: BSD license
* Documentation: https://molssi-seamm.github.io/read_structure_step/index.html
* Code: https://github.com/molssi-seamm/read_structure_step

Features
--------

- Emphasis in ease-of-use and simplicity. The public interface
  consists of a single function.
- Automatic file type recognition.
- OpenBabel as an engine to parse formats, but other engines can be
  easily implemented.
- Easily extensible to new formats.
- Current support for PDB, Mol2 and XYZ files.

Example
+++++++

.. code:: python

 import read_structure_step
 seamm_structure = read_structure_step.read("spc.xyz")

Acknowledgements
----------------

This package was created with Cookiecutter_ and the `molssi-seamm/cookiecutter-seamm-plugin`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`molssi-seamm/cookiecutter-seamm-plugin`: https://github.com/molssi-seamm/cookiecutter-seamm-plugin

Developed by the Molecular Sciences Software Institute (MolSSI_),
which receives funding from the `National Science Foundation`_ under
award ACI-1547580

.. _MolSSI: https://www.molssi.org
.. _`National Science Foundation`: https://www.nsf.gov
