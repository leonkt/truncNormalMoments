## Test environments
- local OSX install, R 4.2.1
- Fedora Linux, R-devel, clang, gfortran (on R-hub)
- (on R-hub)
- win-builder (devel and release)


## R CMD check results

There were no ERRORs or WARNINGs.

There were 3 NOTEs:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Mika Braginsky <mika.br@gmail.com>'

New submission

Possibly misspelled words in DESCRIPTION:
  Fabozzi (12:19)
  Jeffreys (11:43)
  Zhou (11:63)
```

This is a new package. None of these are actually misspellings.


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
