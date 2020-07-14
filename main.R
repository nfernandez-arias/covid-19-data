setwd("~/Insync/nfernand@princeton.edu/Google Drive/PhD - Thesis/Research/covid-19-data")

source("estimation.R")

source("projections.R")

max(UStotals[variable == "projectedDeaths_10"]$value, na.rm = TRUE)
max(UStotals[variable == "deaths"]$value, na.rm = TRUE)

#source("projectionTest.R")
                                