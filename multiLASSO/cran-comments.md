## Test environments
* local macOS install, R 3.6.0
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE:

* checking dependencies in R code ... NOTE
  Unexported object imported by a ':::' call: ‘gimme:::setupConvolve’

  GIMME is an r-package currently maintained internally by the same lab. Updates to the relevant function in GIMME should be used by this package.

## Downstream dependencies
There are currently no downstream dependencies for this package.