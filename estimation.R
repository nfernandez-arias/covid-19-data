# Clear workspace
rm(list = ls())

library(data.table)
library(ggplot2)
library(ggthemr)
ggthemr("flat")# Set theme -- controls all plots
library(gridExtra)

setwd("~/nfernand@princeton.edu/PhD - Thesis/Research/covid-19-data")

#  Parameters for state space estimation of asymptomatic daily R

f0 <- 0.9  # initial fraction asymptomatic
f <- 0.9 # fraction of new cases that are asymptomatic
r <- 0.1 # fraction of cases that resolve each day

worsenFrac <- 0.2 # fraction of asymnptomatic cases that eventually become symptomatic
w <- ((worsenFrac) / (1- worsenFrac)) * r


states <- fread("us-states.csv")

states <- states[fips <= 56]

population <- fread("nst-est2019-alldata.csv")
area <- fread("statesAreas.csv")

setkey(area,State)
setkey(states,state)

states <- area[states]

setnames(states,"State","state")

setkey(population,STATE)
setkey(states,fips)

states <- population[ , .(STATE,POPESTIMATE2019)][states]

setkey(states,state,date)

BigStates <- c("New York","California","Washington","District of Columbia","Illinois","New Jersey","Texas","Michigan","Massachusetts","Florida","Louisiana","Pennsylvania")

states[ , date := as.Date(date)]

states[ , index := rowid(state) - 1]

states[ , logCases := log(cases)]
states[ , logDeaths := log(deaths)]

states[ , newCasesRaw := cases - shift(cases), by = state]
states[ , newDeathsRaw := deaths - shift(deaths), by = state]

states[ , newCases := logCases - shift(logCases), by = state]
states[ , newDeaths := logDeaths - shift(logDeaths), by = state]

states[ , currentCases := cases - deaths]
states[ , newCasesAdjusted := newCasesRaw / shift(currentCases), by = state]

# Con
states[ , newAsymptomaticRaw := newCasesRaw * (f / (1-f)) ]

states[states[ , .I[1], by = state]$V1, newAsymptomaticRaw := (f0 / (1-f0)) * cases]
states[states[ , .I[1], by = state]$V1, newCasesRaw := cases]
states[states[ , .I[1], by = state]$V1, newDeathsRaw := deaths]

filter <- function(cases,newCases,newAsymp,deaths,newDeaths) {
  
  asymp <- vector(mode = "numeric", length = length(cases))
  cases_ongoing <- vector(mode = "numeric", length = length(cases))
  immune <- vector(mode = "numeric", length = length(cases))
  
  asymp[1] <- newAsymp[1]
  
  for (i in 2:length(cases)) {
  
      immune[i] <- immune[i-1] + r * (asymp[i-1] + cases_ongoing[i-1])  
      asymp[i] <- asymp[i-1] * (1 - r - w) + newAsymp[i]
      cases_ongoing[i] <- cases_ongoing[i-1] * (1-r) - newDeaths[i-1] + asymp[i-1] * w + newCases[i]
    
  }
  
  return(list(asymp,cases_ongoing,immune))

}



states[ , c("asymptomatics","cases_ongoing","immune") := filter(cases,newCasesRaw,newAsymptomaticRaw,deaths,newDeathsRaw), by = state]

states[ , vulnerablePopulation := POPESTIMATE2019 - asymptomatics - cases_ongoing - immune - deaths]

states[ , impliedR_asymp := newAsymptomaticRaw / (f * shift(asymptomatics)), by = state]

states[ , impliedT_asymp := newAsymptomaticRaw / (f * shift(vulnerablePopulation) * shift(asymptomatics)), by = state]

fwrite(states,"statesModified.csv")



ggplot(data = states[], aes(x = date, y = log(cases), group = state)) +
  geom_line() +
  facet_wrap(~state) +
  labs(title = "Coronavirus cases by state") + 
  theme(axis.text.x = element_text(angle = 90))


ggsave("US_states_corona_cases.pdf",plot = last_plot(), width = 14, height = 10, units = "in")

ggplot(data = states[], aes(x = date, group = state)) +
  geom_line(aes(y = log(asymptomatics), color = "Asymptomatics")) +
  geom_line(aes(y = log(cases_ongoing), color = "Symptomatics")) + 
  facet_wrap(~state) +
  labs(title = "Coronavirus cases (asymptomatic) by state") + 
  theme(axis.text.x = element_text(angle = 90))

ggsave("US_states_corona_cases_asymptomatic.pdf",plot = last_plot(), width = 14, height = 10, units = "in")

ggplot(data = states[], aes(x = date, y = log(newCasesRaw), group = state)) +
  geom_line() +
  facet_wrap(~state) +
  labs(title = "New Coronavirus cases by state") + 
  theme(axis.text.x = element_text(angle = 90))


ggsave("US_states_corona_cases_new.pdf",plot = last_plot(), width = 14, height = 10, units = "in")

ggplot(data = states[], aes(x = date, y = exp(newCases) -1, group = state)) +
  geom_line() +
  facet_wrap(~state) +
  labs(title = "Implied daily R") + 
  ylim(0,1) + 
  theme(axis.text.x = element_text(angle = 90))


ggsave("US_states_corona_cases_impliedR.pdf",plot = last_plot(), width = 18, height = 13, units = "in")

ggplot(data = states[], aes(x = date, y = newCasesAdjusted, group = state)) +
  geom_line() +
  facet_wrap(~state) +
  labs(title = "Implied daily R") + 
  ylim(0,1) + 
  theme(axis.text.x = element_text(angle = 90))


ggsave("US_states_corona_cases_impliedR_adjusted.pdf",plot = last_plot(), width = 18, height = 13, units = "in")

ggplot(data = states[], aes(x = date, y = impliedR_asymp, group = state)) +
  geom_line() +
  facet_wrap(~state) +
  labs(title = "Implied daily R") + 
  ylim(0,3) + 
  theme(axis.text.x = element_text(angle = 90))


ggsave("US_states_corona_cases_impliedRasymp_adjusted.pdf",plot = last_plot(), width = 18, height = 13, units = "in")

states[ , impliedT_asymp_standardized := (impliedT_asymp - mean(na.omit(impliedT_asymp)) ) / sd(na.omit(impliedT_asymp)), by = state]

ggplot(data = states[], aes(x = date, y = impliedT_asymp * POPESTIMATE2019 / (r + w), group = state)) +
  geom_line() +
  facet_wrap(~state) +
  ylim(0,10) + 
  labs(title = "Number of people infected per infected person (over course of infection)",
       subtitle = "< 1 means past peak (roughly)") + 
  theme(axis.text.x = element_text(angle = 90))

ggsave("US_states_corona_cases_impliedTasymp_adjusted.pdf",plot = last_plot(), width = 18, height = 13, units = "in")


ggplot(data = states[], aes(x = date, y = 1 - 1/(impliedT_asymp * POPESTIMATE2019 / (r + w)), group = state)) +
  geom_line() +
  facet_wrap(~state) +
  ylim(0,1) + 
  labs(title = "Fraction of population infected at peak in ongoing case load",
       subtitle = "holding constant Beta - if value missing, it means the state looked past peak at that point") + 
  theme(axis.text.x = element_text(angle = 90))

ggsave("US_states_corona_cases_fractionInfectedAtPeak.pdf",plot = last_plot(), width = 18, height = 13, units = "in")

ggplot(data = states[], aes(x = date, y = impliedT_asymp * POPESTIMATE2019 / (r + w), group = state)) +
  geom_line() +
  facet_wrap(~state) +
  ylim(0,10) + 
  labs(title = "Implied Beta",
       subtitle = "For each state, this quantity is held constant in projections") + 
  theme(axis.text.x = element_text(angle = 90))


ggsave("US_states_corona_cases_impliedTasymp_adjusted.pdf",plot = last_plot(), width = 18, height = 13, units = "in")




  ggplot(data = states[], aes(x = date, y = log(deaths), group = state)) +
  geom_line() +
  facet_wrap(~state) +
  labs(title = "Coronavirus deaths by state") + 
  theme(axis.text.x = element_text(angle = 90))


ggsave("US_states_corona_deaths.pdf",plot = last_plot(), width = 14, height = 10, units = "in")


ggplot(data = states[state %in% BigStates], aes(x = date, y = log(cases), group = state)) +
  geom_line() +
  facet_wrap(~state) +
  labs(title = "Coronavirus cases by major US state") +
  theme(axis.text.x = element_text(angle = 90))


ggsave("US_bigstates_corona_cases.pdf",plot = last_plot(), width = 11, height = 8, units = "in")


ggplot(data = states[state %in% BigStates], aes(x = date, y = log(deaths), group = state)) +
  geom_line() +
  facet_wrap(~state) +
  labs(title = "Coronavirus deaths by major US state") +
  theme(axis.text.x = element_text(angle = 90))


ggsave("US_bigstates_corona_deaths.pdf",plot = last_plot(), width = 11, height = 8, units = "in")


              


counties <- fread("us-counties.csv")

counties[ , countystate := paste(county,state, sep = ", ")]


BigCounties <- c(c("New York City, New York"),c("San Francisco, California"),
                 c("Los Angeles, California"),
                 c("Cook, Illinois"),
                 c("King, Washington"),
                 c("Mercer, New Jersey"),
                 c("Bergen, New Jersey"),
                 c("Wayne, Texas"),
                 c("Suffolk, Massachusetts"),
                 c("Philadelphia, Pennsylvania"),
                 c("Wayne, Michigan"),
                 c("Miami-Dade, Florida"),
                 c("Orleans, Louisiana"))


counties[ , date := as.Date(date)]

setkey(counties,countystate,date)

counties[ , logCases := log(cases)]
counties[ , logDeaths := log(deaths)]

counties[ , newCasesRaw := cases - shift(cases), by = countystate]
counties[ , newDeathsRaw := deaths - shift(deaths), by = countystate]

counties[ , newCases := logCases - shift(logCases), by = countystate]
counties[ , newDeaths := logDeaths - shift(logDeaths), by = countystate]

counties[ , currentCases := cases - deaths]
counties[ , newCasesAdjusted := newCasesRaw / shift(currentCases), by = countystate]




ggplot(data = counties[countystate %in% BigCounties], aes(x = date, y = log(cases), group = county)) +
  geom_line() +
  facet_wrap(~countystate) +
  labs(title = "Coronavirus cases by major US county",
       subtitle = "Cook = Chicago, King = Seattle, Suffolk = Boston, Wayne = Detroit") + 
  theme(axis.text.x = element_text(angle = 90)) 


ggsave("US_counties_corona_cases.pdf",plot = last_plot(), width = 11, height = 8, units = "in")

ggplot(data = counties[countystate %in% BigCounties], aes(x = date, y = newCasesAdjusted, group = countystate)) +
  geom_line() +
  facet_wrap(~countystate) +
  labs(title = "Implied daily R",
       subtitle = "Cook = Chicago, King = Seattle, Suffolk = Boston, Wayne = Detroit") + 
  ylim(0,1) + 
  theme(axis.text.x = element_text(angle = 90))


ggsave("US_counties_corona_cases_impliedR_adjusted.pdf",plot = last_plot(), width = 18, height = 13, units = "in")


    ggplot(data = counties[countystate %in% BigCounties], aes(x = date, y = log(deaths), group = county)) +
  geom_line() +
  facet_wrap(~countystate) +
  labs(title = "Coronavirus deaths by major US county",
       subtitle = "Cook = Chicago, King = Seattle, Suffolk = Boston, Wayne = Detroit") + 
  theme(axis.text.x = element_text(angle = 90)) 


ggsave("US_counties_corona_deaths.pdf",plot = last_plot(), width = 11, height = 8, units = "in")






