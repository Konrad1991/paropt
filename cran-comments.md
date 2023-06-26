## Resubmission
This is a resubmission. In this version I have:

* Package was archived on CRAN -> Have these issues been solved? 
  * All the issues were solved. 

*The LICENSE file is only needed if you have additional restrictions to the license which you have not? In   that case omit the file and its reference in the DESCRIPTION file. 
  * I removed the LICENSE file and its reference in the DESCRIPTION file.
  
* replace \dontrun{} with \donttest
  * I replaced always \dontrun{} with \donttest{}

## Test environments
* local Ubuntu Linux, R 4.3.0
* win-builder (oldrelease, devel and release)
* for macOS I used devtools::check_mac_release

## R CMD check results
There were no ERRORs or WARNINGs.

There were 2 NOTEs:

* New submission
  Package was archived on CRAN
* the package was previously archived due to a UBSan error. I have since resolved the error.
  - Sanitizer checks using the clang compiler revealed no errors or warnings

* checking installed package size ... NOTE
  installed size is  7.8Mb
  sub-directories of 1Mb or more:
    libs   5.8Mb

## Downstream dependencies

There are currently no downstream dependencies for this package

