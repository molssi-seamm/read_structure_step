=======
History
=======
2024.12.29 -- Bugfix: Issue with reusing systems matching SDF files.
  * The logic was faulty, so if the first structure in an SDF file was a configuration
    of an existing system, it was not added to the system correctly.
    
2024.12.7 -- Minor update due to changes in molsystem
  * The molsystem was updated to improve the handling of properties. This required
    changes in this module to match the new molsystem.
    
2024.11.29 -- Added chemical formula as option for system/configuration names

2024.11.13 -- Bugfix: Issues with the names of properties
  * Fixed syntax of the properties pulled out from MOPAC encoded parameterization
    inputs.
    
2024.11.3 -- Bugfix: MOPAC files with references in comments
  * Fixed a bug that caused a crash when reading MOPAC files with references in the
    comments.
  * Updated the MOPAC reader to the new approach for running MOPAC in the cases that
    it is needed: Z-matrices and mixed inputs that OpenBabel can't handle.
      
2024.8.23 -- Enhancements to directory handling
  * Changed the handling of paths to make them relative to the directory that the step
    is running in. In a loop, for instance, files are relative to the iteration
    directory.
  * When reading files, allow the syntax 'job://<job number>/...' to read files from
    another job. If the job number is omitted, so the prefix is 'job:///...', the top
    directory of the current job is used.

2024.7.28 -- Added new names for systems and configurations
  * Made the naming of systems and configurations consistent with the standard
    parameters for them in the GUI.
  * Removed all "from file" options, replacing them with "title", which means the title
    from the file, if it exists, or the filename if it doesn't.
    
2023.11.16 -- Bugfix: titles in SDF files
  * Crashed reading some SDF files write by SEAMM due to the system and configuration
    names encoded in the title having multiple slashes (/).

2023.11.5 -- Added writers for CIF and mmCIF.

2023.11.2 -- Initial changes to structure handling
  * Moving towards the standard structure handling that has developed across SEAMM.
  
2023.8.30 -- Support for spacegroup symmetry

2023.7.28 -- Implemented ranges for reading XYZ and SDF files.

2023.7.27.1 -- Removed debug printing.

2023.7.27 -- Support for .gz and .bz2 files, and multi-structure .xyz files
  * Handle .gz and .bz2 files for .sdf and .xyz extensions.
  * Handle multi-structure XYZ files with a blank line between records.
    
2023.7.6 -- Bugfixes
  * Fixed output of number of structures written to SDF files, which was off by 1.
  * Cleaned up the output for the write-structure step
    
2023.1.30 -- Fixed issue#43, duplicate systems or configuration created

  * Reading a single structure from e.g. a .sdf file created a second system or
    configuration, depending on the stucture control options.

2023.1.24 -- Added handler for XYZ files and added properties

  * Added a custom handler for XYZ files to cope with some of the variant formats.

    * Supports files with no atom count on the first line

    * Supports the variant used in the Minnesota Solubility database, which has 3 header
      lines.

  * Add capability to store properties into the database for formats such as MOPAC and
    SDF that can handle properties. Also can output the properties when writing the
    files.

  * Fixed bugs if the system name or configuration name is not a string, but rather a number.

2022.10.28 -- Fixed bug reading cif and mmcif files

  * There was a bug that caused a crash when reading cif and mmcif files, and potentially
    some other formats. It has been fixed throughout.

  * The standard error for properties were missing a commma in the property name. The
    comma is standard elsewhere in SEAMM so add it here: '<prop>, stderr'

2022.10.26 -- Handling OpenBabel error messages for MOPAC .mop files
  Hiding messages about errors Kekulizing structures, which doesn't seem to be a serious
  issue, and printing any other messages as warnings.

2021.2.12 (12 February 2021)

  * Updated the README file to give a better description.

  * Updated the short description in setup.py to work with the new installer.

  * Added keywords for better searchability.

2021.2.4 (4 February 2021)
  Updated for compatibility with the new system classes in MolSystem
  2021.2.2 release.

2020.12.5 (5 December 2020)
  Internal: switching CI from TravisCI to GitHub Actions, and in the
  process moving documentation from ReadTheDocs to GitHub Pages where
  it is consolidated with the main SEAMM documentation.

2020.8.1 (1 August 2020)
  Removed leftover debug print statements.

0.9 (15 April 2020)

  * General bug fixing and code cleanup.

  * Part of release of all modules.

0.7.1 (23 November 2019)
  First release on PyPI.
