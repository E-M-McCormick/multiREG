<!-- README.md is generated from README.Rmd. Please edit that file -->
**Network Regularization (netreg)**
===================================

The netreg algorithm is a continually maintained R package.

Program developers are invited to submit changes here at the GitHub
repository.

**The Basics**
==============

-   Missing data is not a problem.

-   Heterogeneous data is not a problem:

    -   No “group” or “common” network map will be forced unless it
        truly describes the majority.

    -   Individual-level nuances will surface after a group or common
        network map is derived (provided one exists).

    -   If desired, subgroups of individuals with similar patterns of
        effects will be generated to aid the researcher in finding
        similar patterns among the varied individual models.

-   Supports a large number of observed variables (capable of estimating
    a whole-brain network map).

    -   Alternatively can support interaction effects between endogenous
        variables, or between exogenous and endogenous variables.

-   Can be freely downloaded by installing the package “netreg” in R.

**Running netreg**
==================

**1. Create two new folders (i.e., directories)**

-   Create a source folder for your time series. This can be anywhere
    that you have permission to read and write.

-   Nothing can be in the source folder other than the time series data.

-   Create an output folder for your results. This must be different
    from the above folder. Alternatively, specifying an output folder
    will cause one to be created if permissions allow.

**2. Extract the time series for your variables**

-   Have each variable be a column, with the rows being the observation
    (e.g., scan in fMRI or a day in daily diary studies).

-   Substitute NA for missing values.

-   Have a separate file for each individual/session.

-   Put each file in the source folder you created in step 1. Do not put
    anything else in this folder.

-   Files must be either comma, space, or tab delimited.

**3. Installing netreg with R**

-   Open an R script and enter into the console:
    `install.packages("netreg")`

-   Once netreg has been installed, you will need to load the package by
    entering: `library(netreg)`

**4. Running gimme**

The *netreg* (or equivalently, *network\_reg*) function requires that
you input:

-   The path to the directory containing your data.

-   How data are separated (e.g., comma-separated values).

-   Whether the data files contain a header row.

All other fields are optional and will go to defaults if no user input
is provided. If no output directory is indicated, all information is
stored as R objects (see tutorial linked above for details).

``` r
output <- netreg(     # can use "netreg" or "network_reg"
  data = '',          # source directory where your data are 
  out = NULL,         # output directory where you'd like your output to go (if NULL, output will only be saved in a list object)
  sep = NULL,         # how data are separated. "" for space; "," for comma, "/t" for tab-delimited
  header = TRUE,      # TRUE or FALSE, is there a header
  ar = TRUE,          # TRUE (default) or FALSE, start with autoregressive paths open
  plot = TRUE,        # TRUE (default) or FALSE, generate plots
  subgroup = FALSE,   # TRUE or FALSE (default), cluster individuals based on similarities in effects
  alpha = .5,         # option to control the elasticnet mixing parameter; alpha = .5 (default), alpha = 1 is the lasso penalty, alpha = 0 is the ridge regression penalty
  penalties = NULL,   # option to specify a matrix of shrinkage parameters that will control the initial search for a group-level network map
  groupcutoff = .75,  # the proportion that is considered the majority at the group level
  subcutoff = .5      # the proportion that is considered the majority at the subgroup level
)        
```

While *netreg* is running you will see information iterate in the
command window. The algorithm will alert you when it is finished.

**Output**
==========

-   The output directory will contain:

    -   **GroupLevel\_PathCountsMatrix**:
