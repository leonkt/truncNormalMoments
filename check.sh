#!/bin/bash

## Usage
##
## Define and pull a Docker image
# export IMAGE="rhub/ubuntu-release"
# export IMAGE="eddelbuettel/r2u:22.04"
# docker pull $IMAGE
##
## Clean up artifact folder
# mkdir _tmp
# rm -rf _tmp/*
##
## Run the check.sh script
# docker run -it --rm --platform linux/amd64 -v $PWD:/home/root/truncnormbayes -w /home/root ${IMAGE:-eddelbuettel/r2u:22.04} /bin/bash truncnormbayes/check.sh
##
## Now check files in the _tmp folder

PKG=truncnormbayes
TMP=$PKG/_tmp
FILE=$TMP/$PKG-check-logs.log
echo $(date) - Script started > $FILE

echo $(date) - Installing package dependencies >> $FILE
apt update -qq && apt install --yes --no-install-recommends pandoc devscripts qpdf
R -q -e "install.packages(c('codetools','roxygen2'))"
R -q -e "remotes::install_deps('$PKG', upgrade='never', dependencies=TRUE)"
# R -q -e "install.packages(c('codetools','roxygen2','rmarkdown','knitr','testthat'))"
# R -q -e "remotes::install_deps('$PKG', upgrade='never', dependencies=NA)"
# R CMD build $PKG
# cp $PKG_*.tar.gz $TMP

echo $(date) - Updating Rd files >> $FILE
R CMD build $PKG
R CMD INSTALL $PKG_*.tar.gz
R -q -e "roxygen2::roxygenize('$PKG', load_code='installed')"
rm $PKG_*.tar.gz

echo $(date) - Building package >> $FILE
R CMD build $PKG
cp $PKG_*.tar.gz $TMP

echo $(date) - Running checks >> $FILE
R CMD check --no-manual --as-cran --output=$TMP $PKG_*.tar.gz
# R_CHECK_DONTTEST_EXAMPLES=false R CMD check --no-manual --as-cran --output=$TMP $PKG_*.tar.gz
# _R_CHECK_DEPENDS_ONLY_=true _R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual --as-cran --output=$TMP $PKG_*.tar.gz

echo $(date) - Done >> $FILE
