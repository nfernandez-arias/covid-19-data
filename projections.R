# Clear workspace
rm(list = ls())

library(data.table)
library(ggplot2)
library(ggthemr)
ggthemr("flat")# Set theme -- controls all plots
library(gridExtra)

setwd("~/nfernand@princeton.edu/PhD - Thesis/Research/covid-19-data")

complete <- tidyr::complete

f0 <- 0.9  # initial fraction asymptomatic    
f <- 0.9 # fraction of new cases that are asymptomatic
r <- 0.05 # fraction of cases that resolve each day

# Factor by which state-specific R is adjusted
factor <- 1

states <- fread("statesModified.csv")


deathRate <- states[ , .SD[.N], by = state][ , sum(na.omit(deaths)) / sum(na.omit(cases))]




states[ , impliedR_asymp_ma3 := (1/5) * Reduce(`+`, shift(impliedR_asymp,n = 0L:4L, type = "lag")), by = state]
states[ , impliedT_asymp_ma3 := (1/5) * Reduce(`+`, shift(impliedT_asymp,n = 0L:4L, type = "lag")), by = state]
Rtable <- states[ , .SD[.N], by = state][ , .(state,impliedR_asymp_ma3)]
Ttable <- states[ , .SD[.N], by = state][ , .(state,impliedT_asymp_ma3)]

setnames(Rtable,"impliedR_asymp_ma3","impliedR_last")
setnames(Ttable,"impliedT_asymp_ma3","impliedT_last")

necessaryR <- r / f

Rtable[ , necessaryFactor:= necessaryR / impliedR_last]

setkey(Rtable,state)
setkey(Ttable,state)
d = (deathRate / (1 - deathRate)) * r # Projected death rate per day  

states[ , date := as.Date(date)]

# Extend dataset for each state until end of projection period
temp <- states[ , .(date = seq.Date(min(date),as.Date("2020-08-01"), by = "day")), by = state]

setkey(temp,state,date)
setkey(states,state,date)
  
states <- states[temp]

states <- Rtable[states]
states <- Ttable[states]


# Make projections

project <- function(asymptomatics,cases,deaths,vulnerablePopulation,impliedT)  {

  projectedNewAsymptomatics <- vector(mode = "numeric", length = length(asymptomatics))
  projectedNewCases <- vector(mode = "numeric", length = length(asymptomatics))
  projectedNewDeaths <- vector(mode = "numeric", length = length(asymptomatics))
  
  projectedAsymptomatics <- vector(mode = "numeric", length = length(asymptomatics))
  projectedCases <- vector(mode = "numeric", length = length(asymptomatics))
  projectedDeaths <- vector(mode = "numeric", length = length(asymptomatics))
  projectedVulnerablePopulation <- vector(mode = "numeric", length = length(asymptomatics))
  
  iMin <- max(which(!is.na(cases))) + 1
  
  projectedAsymptomatics[iMin - 1] <- asymptomatics[iMin - 1]
  projectedCases[iMin - 1] <- cases[iMin - 1]
  projectedDeaths[iMin -1] <- deaths[iMin - 1]
  projectedVulnerablePopulation[iMin -1] <- vulnerablePopulation[iMin -1]
  
  for (i in iMin:length(asymptomatics)) {
    
    projectedNewAsymptomatics[i] <- f * factor * impliedT[i] * projectedVulnerablePopulation[i-1] * projectedAsymptomatics[i-1]
    projectedAsymptomatics[i] <- projectedAsymptomatics[i-1] * (1-r) + projectedNewAsymptomatics[i]
    
    projectedNewCases[i] <- (1-f) * factor * impliedT[i] * projectedVulnerablePopulation[i-1] * projectedAsymptomatics[i-1]
    projectedCases[i] <- projectedCases[i-1] * (1-r) + projectedNewCases[i]
    
    projectedNewDeaths[i] <- d * projectedCases[i]
    projectedDeaths[i] <- projectedDeaths[i-1] + projectedNewDeaths[i]

    projectedVulnerablePopulation[i] <- projectedVulnerablePopulation[i-1] - projectedNewAsymptomatics[i] - projectedNewCases[i]
        
  }
  
  return(list(projectedNewAsymptomatics,projectedAsymptomatics,projectedNewCases,projectedCases,projectedNewDeaths,projectedDeaths,projectedVulnerablePopulation))
  
}

states[ , c("projectedNewAsymptomatics","projectedAsymptomatics","projectedNewCases","projectedCases","projectedNewDeaths","projectedDeaths","projectedVulnerablePopulation") := project(asymptomatics,cases_ongoing,deaths,vulnerablePopulation,impliedT_last), by = state ]




# Make plots

states[ projectedCases == 0, projectedCases := NA]
states[ projectedAsymptomatics == 0, projectedAsymptomatics := NA]
states[ projectedDeaths == 0, projectedDeaths := NA]

states[ cases_ongoing == 0, cases_ongoing := NA]
states[ asymptomatics == 0, asymptomatics := NA]
#states[ projectedDeaths == 0.0, deaths := NA]

UStotals <- states[ , .(projectedCases = sum(na.omit(projectedCases)),projectedAsymptomatics = sum(na.omit(projectedAsymptomatics)), projectedDeaths = sum(na.omit(projectedDeaths)),
                        cases_ongoing = sum(na.omit(cases_ongoing)), asymptomatics = sum(na.omit(asymptomatics)), deaths = sum(na.omit(deaths))), by = date]

setkey(UStotals,date)

UStotals[ projectedCases == 0, projectedCases := NA]
UStotals[ projectedAsymptomatics == 0, projectedAsymptomatics := NA]
UStotals[ projectedDeaths == 0, projectedDeaths := NA]

UStotals[ cases_ongoing == 0, cases_ongoing := NA]
UStotals[ asymptomatics == 0, asymptomatics := NA]
UStotals[ deaths == 0, deaths := NA]

ggplot(UStotals, aes(x = date)) + 
  geom_line(aes(y = log(cases_ongoing), linetype = "Actual", color = "Cases (ongoing)")) + 
  geom_line(aes(y = log(projectedCases), linetype = "Projected", color = "Cases (ongoing)")) +
  geom_line(aes(y = log(asymptomatics), linetype = "Actual", color = "Asymptomatics (ongoing)")) + 
  geom_line(aes(y = log(projectedAsymptomatics), linetype = "Projected", color = "Asymptomatics (ongoing)")) +
  geom_line(aes(y = log(deaths), linetype = "Actual", color = "Deaths (cumulative)")) + 
  geom_line(aes(y = log(projectedDeaths), linetype = "Projected", color = "Deaths (cumulative)")) + 
  scale_x_date(date_breaks = "2 weeks") + 
  labs(title = "US total") + 
  ylab("Log(# of people)") + 
  xlab("Date") + 
  theme(axis.text.x = element_text(angle = 90))

ggsave("US_corona_projections.pdf",plot = last_plot(), width = 10, height = 6, units = "in")


ggplot(states, aes(x = date)) + 
  geom_line(aes(y = log(cases_ongoing), linetype = "Ongoing cases (actual, estimated)")) + 
  geom_line(aes(y = log(projectedCases), linetype = "Ongoing cases (projected)")) + 
  facet_wrap(~state) + 
  theme(axis.text.x = element_text(angle = 90))

ggsave("US_states_corona_cases_projections.pdf",plot = last_plot(), width = 18, height = 13, units = "in")

ggplot(states, aes(x = date)) + 
  geom_line(aes(y = log(asymptomatics), linetype = "Asymptomatics (actual, estimated)")) + 
  geom_line(aes(y = log(projectedAsymptomatics), linetype = "Asymptomatics (projected)")) + 
  facet_wrap(~state) + 
  theme(axis.text.x = element_text(angle = 90))

ggsave("US_states_corona_asymptomatics_projections.pdf",plot = last_plot(), width = 18, height = 13, units = "in")


ggplot(states, aes(x = date)) + 
  geom_line(aes(y = log(deaths), linetype = "Deaths (actual)")) + 
  geom_line(aes(y = log(projectedDeaths), linetype = "Deaths (projected)")) + 
  facet_wrap(~state) + 
  theme(axis.text.x = element_text(angle = 90))

ggsave("US_states_corona_deaths_projections.pdf",plot = last_plot(), width = 18, height = 13, units = "in")


# Plot necessary factors
ggplot(Rtable[ , .(necessaryFactor)] , aes(x = necessaryFactor)) + 
  geom_histogram(bins = 50) + 
  stat_density(geom = "line")

ggsave("necessaryR_adjustment.pdf", plot = last_plot(), width = 10, height = 6, units = "in")





