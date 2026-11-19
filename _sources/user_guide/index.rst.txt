.. _user-guide:

**********
User Guide
**********
The Read Structure plug-in reads structures from a wide range of molecular file formats.

Files with several structures
=============================
Many formats, such as SDF and extended XYZ files, can hold more than
one structure. By default the step reads all of them. **Structures to read**
(``indices``) chooses which: a comma-separated list of indices and ranges
``start:stop[:step]``, counting from 1 with the stop included, e.g. ``1`` for only the
first structure, or ``1:10:2, 20:end``.

The first structure read is handled as **Structure handling** says: by default it
overwrites the current configuration. Each later structure is handled as **Subsequent
structure handling** says: by default it goes into a new system and configuration.

The last structure read becomes the current configuration, and the steps that follow
work on it. When the step reads more than one structure, its output says which one is
now current, e.g.::

    The current configuration, which the following steps use, is now structure 2
    of 2, C4 H8 O:
               System name = 109-92-2_0
        Configuration name = C=COCC

Check this line when a file may hold more structures than you expect. If only the
first is wanted, set **Structures to read** to ``1``.

Plain XYZ files with several structures
---------------------------------------
Plain XYZ files come in several layouts, and the step reads each of them, one after
another in the same file:

* The standard one: the number of atoms on the first line, a comment (often the name)
  on the second, then one line per atom. The structures may follow one another
  directly or be separated by blank lines.
* Without the number of atoms: a blank line in its place, or no line at all, then the
  comment and the atoms. Here a blank line separates the structures.
* The Minnesota Solvation Database's layout: the comment, a blank line, the charge and
  spin multiplicity, then the atoms, with a blank line between structures. The charge
  and multiplicity are kept.

Set **System name** and **Configuration name** to ``title`` to name each structure
after its comment line.

MOPAC files
===========
A MOPAC input file (``.mop``) is read as MOPAC itself reads its geometry: Cartesian
coordinates, internal coordinates (a Z-matrix, with dummy atoms), or a mix of the two.
MOPAC is not needed. The charge and spin multiplicity come from the keywords
(``CHARGE=``, ``SINGLET`` ... ``NONET``, ``MS=``), and the title line can name the
system and configuration.

The bonds are perceived from the distances, and their orders and the formal charges are
assigned molecule by molecule. Metal complexes are ionic + dative: every metal-ligand
bond is a dative bond from the ligand's donor atom to the metal, each ligand has its
closed-shell Lewis structure (CO neutral, Cp-, Cl-, acetate), metal-metal bonds are
single bonds, and the metals carry the rest of the charge: Fe(CO)5 is Fe(0), Cp2TiCl2
Ti4+. The alkali and alkaline-earth metals are left unbonded, as ions. Distance-based
perception can be wrong for long metal-ligand contacts, agostic hydrogens, ion pairs and
clusters, so check such structures. A geometry with translation vectors (``Tv``) is
read with Open Babel, as before, and needs MOPAC when it has dummy atoms.

Reference data
--------------
The reference data of a MOPAC parameterization file -- the values on its reference
line and their sources -- are stored as properties of the configuration, in SEAMM's
units, with the names of the PM-8 reference set's SDF records:

* ``H=``, ``S=``, ``CP=``, ``D=`` and ``I=`` (or ``IA=``, ``IE=``, ``IP=``, ``I1=``) as
  ``enthalpy of formation``, ``entropy``, ``constant pressure heat capacity``,
  ``dipole moment`` and ``ionization energy``, in kJ/mol, J/mol/K, D and eV, followed
  by ``#experiment`` -- or ``#<source>`` for a calculated source such as PW91D, or
  ``#estimate``;
* for each, ``,uncertainty`` (from ``H=value,uncertainty``), ``,reference`` (the
  source's key, from ``HR=``, ``DR=``, ``IR=``), ``,weight`` (``HWT=`` ...) and
  ``,relative to`` (the species of a relative value), e.g. ``enthalpy of
  formation,reference#experiment``;
* ``geometry,reference``, ``geometry,labels`` and ``geometry,weight`` for a geometric
  reference (``GR=``);
* ``enthalpy of formation#MOPAC#PM7`` for the heat of formation in the comment lines;
* ``MOPAC keywords``, ``MOPAC reference line``, ``MOPAC file``, and ``MOPAC geometry``
  when Cartesian coordinates would not keep everything the geometry block holds, so
  that the file can be written again as it was. Its blank lines are written as ``.``,
  since a data item of an SDF file cannot hold an empty line.

The sources' keys refer to a BibTeX bibliography, ``references.bib``. If there is one in
the file's directory or a parent directory, each value also gets ``,locator`` (where in
the source it is) and ``,note`` (how it was derived) from the entry's
``mopac_locator`` and ``mopac_derivation`` fields; the cited entries are kept, as
written, in the ``bibliography`` property; and a source whose entry has
``mopac_kind = {calculated}`` tags its values with its key rather than ``experiment``
(``{estimate}``: ``estimate``). Without a bibliography every source counts as
experimental.

..
   The following sections cover accessing and controlling this functionality.

   .. toctree::
      :maxdepth: 2
      :titlesonly:

Index
=====

* :ref:`genindex`
