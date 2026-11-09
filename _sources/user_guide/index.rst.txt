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

..
   The following sections cover accessing and controlling this functionality.

   .. toctree::
      :maxdepth: 2
      :titlesonly:

Index
=====

* :ref:`genindex`
