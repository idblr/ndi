## This is the sixth resubmission

* Actions taken since previous submission:
  * 'DescTools' is now Suggests to fix Rd cross-references NOTE
  * Fixed 'lost braces in \itemize' NOTE for `anthopolos()`, `atkinson()`, `bell()`, `bemanian_beyer()`, `bravo()`, `duncan()`, `krieger()`, `messer()`, `powell_wiley()`, `sudano()`, and `white()` functions
  * Fixed 'Moved Permanently' content by replacing the old URL with the new URL
  * Fixed citation for Slotman _et al._ (2022) in CITATION

* Documentation for DESCRIPTION, README, NEWS, and vignette references the following DOIs, which throws a NOTE but are a valid URL:
  * <https://doi.org/10.1111/j.1749-6632.2009.05333.x>
  * <https://doi.org/10.2307/2223319>
  * <https://doi.org/10.2307/2088328>
  * <https://doi.org/10.2307/270845>
  * <https://doi.org/10.1080/17445647.2020.1750066>
  * <https://doi.org/10.2307/3644339>
  * <https://doi.org/10.2307/2084686>
  
* Some tests and examples for `anthopolos()`, `atkinson()`, `bell()`, `bemanian_beyer()`, `bravo()`, `duncan()`, `gini()`, `krieger()`, `messer()`, `powell_wiley()`, `sudano()`, and `white()` functions require a Census API key so they are skipped if NULL or not run

## Test environments
* local Windows install, R 4.2.1
* win-builder, (devel, release, oldrelease)
* Rhub
  * Fedora Linux, R-devel, clang, gfortran
  * Ubuntu Linux 20.04.1 LTS, R-release, GCC
  * Windows Server 2022, R-devel, 64 bit
  * Windows Server 2008 R2 SP1, R-release, 32⁄64 bit
  * Oracle Solaris 10, x86, 32 bit, R-release
  * macOS 10.13.6 High Sierra, R-release, CRAN's setup

## R CMD check results
0 errors | 0 warnings | 0 notes

## Submitted by Maintainer
