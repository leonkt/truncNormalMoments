## Test environments

The package was tested with `R CMD check --as-cran` on the following platforms:

* Ubuntu Linux GCC (oldrel, release, devel),
* Windows x86_64  (oldrel, release, devel),
* MacOS (release).

## R CMD check results

There were no ERRORs or WARNINGs. The following NOTEs were displayed:

```
* checking CRAN incoming feasibility ... [6s/15s] NOTE
Maintainer: ‘Peter Solymos <peter@analythium.io>’

New maintainer:
  Peter Solymos <peter@analythium.io>
Old maintainer(s):
  Mika Braginsky <mika.br@gmail.com>
```


The package has new maintainer: Peter Solymos.

```
* checking installed package size ... NOTE
  installed size is 41.2Mb
  sub-directories of 1Mb or more:
    libs  41.0Mb
```

The large installed size is toe to the rstan compiled code being included.

```
* checking for GNU extensions in Makefiles ... NOTE
GNU make is a SystemRequirements.
```

GNU make is required by rstan.
