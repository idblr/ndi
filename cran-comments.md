## This is the first resubmission

* Actions taken since previous submission based on feedback from Beni Altmann:
  * Fixed invalid URL and typos in package README.md

* Documentation for DESCRIPTION and README references the following DOI, which throws a NOTE but is a valid URL:
  * <https://doi.org/10.1111/j.1749-6632.2009.05333.x>
  
* Some tests and examples for `messer()` and `powell_wiley()` functions require require a Census API key so they are skipped if NULL or not run

## Test environments
* local OS X install, R 4.2.1
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
