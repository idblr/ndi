## This is the eighth resubmission

* Actions taken since previous submission:
  * Fixed broken URLs in 'theil.Rd', 'ndi-package.Rd', 'ndi2.html', README, and NEWS

* Words that throw a NOTE by DEBIAN and WINDOWS pre-tests as possibly misspelled but are OK: "geospatial" 
* The win-builder oldrelease throws a NOTE that "Author field differs from that derived from Authors@R". The behavior is OK because ORCID has different formatting but same information
  
* Some tests and examples for `anthopolos()`, `atkinson()`, `bell()`, `bemanian_beyer()`, `bravo()`, `denton()`, `denton_cuzzort()`, `duncan()`, `duncan_cuzzort()`, `duncan_duncan()`, `gini()`, `hoover()`, `james_taeuber()`, `krieger()`, `lieberson()`, `massey()`, `massey_duncan()`, `messer()`, `powell_wiley()`, `sudano()`, `theil()`, `white()`, and `white_blau()` functions require a Census API key so they are skipped if NULL or not run

## Test environments
* local Windows install, R 4.5.1
* win-builder, (devel, release, oldrelease)
* R-CMD-check on GitHub
  * macos-latest (release)
  * windows-latest (release)
  * ubuntu-latest (devel)
  * ubuntu-latest (release)
  * ubuntu-latest (oldrel-1)
* Rhub v2
  * macos-15 on GitHub, ASAN + UBSAN on macOS (`m1-san`)
  * macos-13 on GitHub(`macos`)
  * Fedora Linux 40 (Container Image) (`gcc-asan`)
  * Ubuntu 22.04.5 LTS (`ubuntu-clang`)
  * Ubuntu 22.04.5 LTS (`ubuntu-gcc12`)

## R CMD check results
0 errors | 0 warnings | 0 notes

## Submitted by Maintainer
