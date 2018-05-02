## calculate_scores.R

## This script calculates OHI scores with the `ohicore` package.
## - configure_toolbox.r ensures your files are properly configured for `ohicore`.
## - The `ohicore` function CalculateAll() calculates OHI scores.

## set working directory for all OHI calculations
setwd("~/github/toolbox-demo/region2017")

## run the configure_toolbox.r script to check configuration
source("configure_toolbox.R")

## calculate scenario scores
scores <- ohicore::CalculateAll(conf, layers)

## save scores as scores.csv
readr::write_csv(scores, 'scores.csv', na='')
