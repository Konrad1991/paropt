## Test environments
* local Ubuntu Linux, R 4.2.2
* win-builder (oldrelease, devel and release)
* rhub: "macos-highsierra-release"      "macos-highsierra-release-cran"

## R CMD check results
There were no ERRORs or NOTEs.

1 Warning for:

D:/RCompile/CRANpkg/lib/4.2/ast2ast/include/etr_bits/vec.hpp:46:7: warning: 'etr::VEC<double>::nrows' will be initialized after [-Wreorder]
D:/RCompile/CRANpkg/lib/4.2/ast2ast/include/etr_bits/vec.hpp:45:7: warning: 'int etr::VEC<double>::ncols' [-Wreorder]
D:/RCompile/CRANpkg/lib/4.2/ast2ast/include/etr_bits/vec.hpp:101:3: warning:   when initialized here [-Wreorder]

This is due to code imported from the package 'ast2ast'. As ncols and nrows are independently initialized the order is not important here.

## Downstream dependencies

There are currently no downstream dependencies for this package

