=======
History
=======

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
