# Clear workspace
rm(list = ls())

library(data.table)
library(ggplot2)
library(ggthemr)
ggthemr("flat")# Set theme -- controls all plots
library(gridExtra)

setwd("~/nfernand@princeton.edu/PhD - Thesis/Research/covid-19-data")

states <- fread("us-states.csv")

setkey(states,state,date)

BigStates <- c("New York","California","Washington","District of Columbia","Illinois","New Jersey","Texas","Michigan","Massachusetts","Florida","Louisiana","Pennsylvania")

states[ , date := as.Date(date)]

states[ , logCases := log(cases)]
states[ , logDeaths := log(deaths)]

states[ , newCasesRaw := cases - shift(cases), by = state]
states[ , newDeathsRaw := deaths - shift(deaths), by = state]

states[ , newCases := logCases - shift(logCases), by = state]
states[ , newDeaths := logDeaths - shift(logDeaths), by = state]

states[ , currentCases := cases - deaths]
states[ , newCasesAdjusted := newCasesRaw / shift(currentCases), by = state]

ggplot(data = states[], aes(x = date, y = log(cases), group = state)) +
  geom_line() +
  facet_wrap(~state) +
  labs(title = "Coronavirus cases by state") + 
  theme(axis.text.x = element_text(angle = 90))


ggsave("US_states_corona_cases.pdf",plot = last_plot(), width = 14, height = 10, units = "in")

ggplot(data = states[], aes(x = date, y = log(newCasesRaw), group = state)) +
  geom_line() +
  facet_wrap(~state) +
  labs(title = "New Coronavirus cases by state") + 
  theme(axis.text.x = element_text(angle = 90))


ggsave("US_states_corona_cases.pdf",plot = last_plot(), width = 14, height = 10, units = "in")

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

counties[ , newCases := cases - shift(cases, n = 1L, type = "lag"), by = state]
counties[ , newDeaths := deaths - shift(deaths, n = 1L, type = "lag"), by = state]

ggplot(data = counties[countystate %in% BigCounties], aes(x = date, y = log(cases), group = county)) +
  geom_line() +
  facet_wrap(~countystate) +
  labs(title = "Coronavirus cases by major US county",
       subtitle = "Cook = Chicago, King = Seattle, Suffolk = Boston, Wayne = Detroit") + 
  theme(axis.text.x = element_text(angle = 90)) 


ggsave("US_counties_corona_cases.pdf",plot = last_plot(), width = 11, height = 8, units = "in")


ggplot(data = counties[countystate %in% BigCounties], aes(x = date, y = log(deaths), group = county)) +
  geom_line() +
  facet_wrap(~countystate) +
  labs(title = "Coronavirus deaths by major US county",
       subtitle = "Cook = Chicago, King = Seattle, Suffolk = Boston, Wayne = Detroit") + 
  theme(axis.text.x = element_text(angle = 90)) 


ggsave("US_counties_corona_deaths.pdf",plot = last_plot(), width = 11, height = 8, units = "in")

