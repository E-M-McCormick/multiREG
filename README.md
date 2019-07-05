<!-- README.md is generated from README.Rmd. Please edit that file -->

# **multiLASSO: Directed Connectivity Search Using Regularized Regression**

The multiLASSO algorithm is a continually maintained R package.

Program developers are invited to submit changes here at the GitHub
repository.

# **The Basics**

  - Missing data is not a problem.

  - Heterogeneous data is not a problem:
    
      - No “group” or “common” network map will be forced unless it
        truly describes the majority.
    
      - Individual-level nuances will surface after a group or common
        network map is derived (provided one exists).
    
      - If desired, subgroups of individuals with similar patterns of
        effects will be generated to aid the researcher in finding
        similar patterns among the varied individual models.

  - Supports a large number of observed variables (capable of estimating
    a whole-brain network map).
    
      - **NOTE:** Large numbers of variables will significantly increase
        run time.
    
      - Alternatively can support interaction effects between endogenous
        variables, or between exogenous and endogenous variables.

  - Can be freely downloaded by installing the package “multiLASSO” in
    R.

# **Running multiLASSO**

**1. Create two new folders (i.e., directories)**

  - Create a source folder for your time series. This can be anywhere
    that you have permission to read and write.

  - Nothing can be in the source folder other than the time series data.

  - Create an output folder for your results. This must be different
    from the above folder. Alternatively, specifying an output folder
    will cause one to be created if permissions allow.

**2. Extract the time series for your variables**

  - Have each variable be a column, with the rows being the observation
    (e.g., scan in fMRI or a day in daily diary studies).

  - Substitute NA for missing values.

  - Have a separate file for each individual/session.

  - Put each file in the source folder you created in step 1. Do not put
    anything else in this folder.

  - Files must be either comma, space, or tab delimited.

**3. Installing multiLASSO with R**

  - Open an R script and enter into the console:
    `install.packages("multiLASSO")`

  - Once multiLASSO has been installed, you will need to load the
    package by entering: `library(multiLASSO)`

**4. Running multiLASSO**

The *multiLASSO* function requires that you input:

  - The path to the directory containing your data.

  - How data are separated (e.g., comma-separated values).

  - Whether the data files contain a header row.

All other fields are optional and will go to defaults if no user input
is provided. If no output directory is indicated, all information is
stored as R objects (see tutorial linked above for details).

``` r
output <- multiLASSO(
  data = '',          # source directory where your data are 
  out = NULL,         # output directory where you'd like your output to go (if NULL, output will only be saved in a list object)
  sep = NULL,         # how data are separated. "" for space; "," for comma, "/t" for tab-delimited
  header = TRUE,      # TRUE or FALSE, is there a header
  ar = TRUE,          # TRUE (default) or FALSE, start with autoregressive paths open
  plot = TRUE,        # TRUE (default) or FALSE, generate plots
  alpha = .5,         # option to control the elasticnet mixing parameter; alpha = .5 (default), alpha = 1 is the lasso penalty, alpha = 0 is the ridge regression penalty
  model_crit = 'bic', # model critera to use in model selection (default to BIC); other options include 'cv' (cross-validation), 'aic', 'aicc', 'hqc'
  penalties = NULL,   # option to specify a matrix of shrinkage parameters that will control the initial search for a group-level network map
  groupcutoff = .75,  # the proportion that is considered the majority at the group level
)        
```

While *multiLASSO* is running you will see information iterate in the
command window. The algorithm will alert you when it is finished.

# **Output**

  - The output directory will contain:
    
      - **function\_summary**: Contains details of function arguments
        and variable names (grouped by variable type) used in the
        current run.
    
      - **indivPathEstimates**: Contains details (iv & dv), regression
        estimate, and level (e.g., group, individual) for each path and
        each individual.
    
      - **groupPathCountsMatrix**: Contains counts of total number of
        paths, including contemporaneous, lagged, and interactions,
        estimated for the sample. The row variable is the predictor and
        the column variable is the outcome variable.
    
      - **groupPathCountsPresent**: Contains matrix of group path
        identities (1 = group-level path, 0 = other). The row variable
        is the predictor and the column variable is the outcome
        variable.
    
      - **groupMainEffectsPlots**: Produced if plot = TRUE. Contains
        figure with group, subgroup (if subgroup = TRUE), and
        individual-level paths for the sample. Black paths are
        group-level, green paths are subgroup-level, and grey paths are
        individual-level, where the thickness of the line represents the
        count. Contemporaneous paths are solid and lagged paths are
        dashed.
    
      - **groupInteractionsPlots**: Produced if plot = TRUE. Contains
        figure with group, subgroup (if subgroup = TRUE), and
        individual-level interaction paths for the sample. Black paths
        are group-level, green paths are subgroup-level, and grey paths
        are individual-level, where the thickness of the line represents
        the count. Interactions are represented by a knot combining the
        relevant predictors, with an arrow from the knot indicating the
        relevant outcome variable. Will only be created if interactions
        are specified.

  - In individual output directory (*where id represents the original
    file name for each individual*):
    
      - ***id*Betas**: Contains individual-level estimates of each path
        for each individual. Rows indicte predictors, columns indicate
        outcomes.
    
      - ***id*MainEffectsPlot**: Contains individual-level plots for
        main effects paths. Red paths represent positive weights and
        blue paths represent negative weights. Contemporaneous paths are
        solid and lagged paths are dashed.
    
      - ***id*InteractionsPlot**: Contains individual-level plots for
        interaction paths. Interactions are represented by a knot
        combining the relevant predictors, with an arrow from the knot
        indicating the relevant outcome variable. Red paths represent
        positive weights and blue paths represent negative weights.
        Contemporaneous paths are solid and lagged paths are dashed.
        Will only be created if interactions are specified.

# **FAQ**

**How many time points do I need?** This is a difficult question since
it will be related to the number of variables you are using. Rules of
thumb for any analysis can generally be used: the more the better\!
Having at least 100 time points is recommended, but adequate results
have been obtained in simulation studies with only T = 60.

**Do all individuals have to have the same number of observations (T)?**
No.

**How many people do I need in my sample?** For *multiLASSO*, reliable
results are obtained with as few as 10 participants. Remember that in
this context, power to detect effects is determined by the number of
time points rather than the number of individuals. Still, having at
least 10 individuals helps *multiLASSO* to detect signal from noise by
looking for effects that consistently occur.

**What do I do if I obtain an error?** Do some initial trouble-shooting.

1.  Ensure that all of your individuals have the same number of
    variables (columns) in their data sets.

2.  Ensure that all variables have variability (i.e., are not constant).
    *multiLASSO* will let you know if this is the case.

3.  Ensure your path directories are correct.

4.  Ensure that the columns are variables and the rows contain the
    observations across time.

5.  If all of this is correct, please email the error you received, code
    used to run *multiLASSO*, and the data (we promise not to use it or
    share it) to: <gimme@unc.edu>.
