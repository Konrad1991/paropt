This release removes the linking to libatomic. Which lead to problems on the Windows/aarch64 platform.

## Test environments
* local Ubuntu Linux, R 4.3.2
* win-builder (oldrelease, devel and release)
* for macOS I used devtools::check_mac_release

## R CMD check results
There were no ERRORs or WARNINGs.

There were 1 NOTEs:

* checking installed package size ... NOTE
  installed size is  8.0Mb
  sub-directories of 1Mb or more:
    libs   6.0Mb

## Downstream dependencies

There are currently no downstream dependencies for this package

