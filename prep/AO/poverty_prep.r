## poverty_prep.r

## In this example we will prepare local data as a new data layer `ao_poverty`. The 'local data' will be dummy data in this example. We will import data that is already formatted as a data layer for the Toolbox, save it in the "layers" folder, and finally make sure it is properly registered in "layers.csv".

## ***************************************************************************** ##

## setup: libraries, file paths ----
library(tidyverse) # install.packages('tidyverse')
dir_layers <- file.path('~/github/toolbox-demo/region2016/layers')


## import dummy 'local data' that is already formatted nicely. Note the naming convention of the data file: it is "goalcode_layername_assessmentYEAR.csv".
data_file  <- file.path(dir_layers, 'ao_need_gl2016.csv')
d <- readr::read_csv(data_file)

## look at a summary to see the range of values
summary(d)


## make dummy changes to data by adding random numbers within a small range. This will keep the values between 0-1, which is probably what the models are expecting, and we won't be changing them at this point.
d2 <- d %>%
  mutate(value = runif(value, min=0, max=0.2))


## save this local data layer in "layers" folder with the same naming convention as above format

readr::write_csv(d2, file.path(dir_layers, "ao_poverty_demo2017.csv"))


## You will see a notification in the "Git" window in RStudio once the layer is saved successfully. Now we need to register this layer in "layers.csv". Let's go back to the tutorial and see how to do that.
