## install_ohicore.r

## install_ohicore.r installs the most current version of the `ohicore` R package.
## You will need to only do this once. `ohicore` requires recent versions of R and RStudio.
  ## 1. Update R. Download the latest version at http://cran.r-project.org.
  ## 2. Update RStudio. Download the latest version at http://www.rstudio.com/products/rstudio/download
  ## 3. Source or run this script.

## delete any existing version of `ohicore`
for (p in c('ohicore')){
  if (p %in% rownames(installed.packages())){
    lib = subset(as.data.frame(installed.packages()), Package==p, LibPath, drop=TRUE)
    remove.packages(p, lib)
  }
}

## install dependencies
for (p in c('devtools', 'git2r')){
  if (!require(p, character.only=TRUE)){
    install.packages(p); require(p, character.only=T)
  }
}

## install most current version of ohicore from github. Warnings are fine, but make sure there are no errors.
devtools::install_github('ohi-science/ohicore')

