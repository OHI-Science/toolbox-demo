## load any libraries needed across website pages
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(knitr)
})

## brewed vars
study_area      <- "Toolbox Demo"
key             <- "toolbox-demo"
dir_scenario_gh <- "https://raw.githubusercontent.com/OHI-Science/toolbox-demo/master/region2016"

## knitr options for all webpages
knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE)

## read in variables if they exist (i.e. don't try for prep repos)
scores_csv <- file.path(dir_scenario_gh, 'scores.csv')
layers_csv <- file.path(dir_scenario_gh, 'layers.csv')
conf_csv   <- file.path(dir_scenario_gh, 'conf/goals.csv')

if (RCurl::url.exists(scores_csv)) scores <- readr::read_csv(scores_csv)
if (RCurl::url.exists(layers_csv)) layers <- readr::read_csv(layers_csv)
if (RCurl::url.exists(conf_csv))   weight <- readr::read_csv(conf_csv) %>%
  select(goal, weight)

## save local copy of conf/goals.Rmd
conf_goals_rmd <- file.path(dir_scenario_gh, 'conf/web/goals.Rmd')

if (RCurl::url.exists(scores_csv)) {
  conf_goals <- readr::read_lines(conf_goals_rmd)
  readr::write_lines(conf_goals, path = 'conf_goals.Rmd', append = FALSE)
}

