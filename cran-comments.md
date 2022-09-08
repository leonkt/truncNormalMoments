## Test environments
- macOS 12.2, R 4.2.1 (local install)
- Ubuntu 20.04, R 4.2.1 (on R-hub)
- win-builder (devel and release)


## R CMD check results

There were no ERRORs or WARNINGs.

There were 2 NOTEs:


```
* checking dependencies in R code ... NOTE
Namespaces in Imports field not imported from:
  'RcppParallel' 'rstantools'
  All declared Imports should be used.
```

These packages are all used.


```
* checking for GNU extensions in Makefiles ... NOTE
GNU make is a SystemRequirements.
```

GNU make is required by rstan.
