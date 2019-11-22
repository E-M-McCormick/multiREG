## Test environments
* local macOS install, R 3.6.0
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs.

## Downstream dependencies
There are currently no downstream dependencies for this package.

## Resubmission
This is a resubmission.

In this version, I have included @param documentation for all functions, included an executable example in the main function, removed setwd from all functions, described all returned objects. I have also changed the name to avoid confusion with other methods.

Additionally, I have removed all commands to change the r directory, included a verbose option to the main function which suppresses all print functions if desired, explained all acronymns in description text, and unwraped the toy example in the documentation. I have also included a link to documentation on the multiREG github site.

I have explicated the GIMME acroynm and provided a reference for further information.