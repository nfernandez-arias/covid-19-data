# Clear workspace
rm(list = ls())

library(data.table)
library(ggplot2)
library(ggthemr)
ggthemr("flat")# Set theme -- controls all plots
library(gridExtra)

setwd("~/nfernand@princeton.edu/PhD - Thesis/Research/covid-19-data")

complete <- tidyr::complete
na.locf <- zoo::na.locf

f0 <- 0.9  # initial fraction asymptomatic    
f <- 0.9 # fraction of new cases that are asymptomatic
r <- 0.1 # fraction of cases that resolve each day

worsenFrac <- 0.2 # fraction of asymnptomatic cases that become symptomatic
w <- ((worsenFrac) / (1- worsenFrac)) * r


# Factor by which state-specific R is adjusted
factor <- 1

states <- fread("statesModified.csv")


deathRateStates <- states[ , .SD[.N], by = state][ , .(deathRate = sum(na.omit(deaths)) / sum(na.omit(cases))), by = state]
deathRate <- states[ , .SD[.N], by = state][ , sum(na.omit(deaths)) / sum(na.omit(cases))]

setkey(deathRateStates,state)
setkey(states,state)

states <- deathRateStates[states]

states[ , deathRate := (0.5 * deathRate / (1- 0.5 * deathRate)) * r]

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
d =   (deathRate / (1 - deathRate)) * r # Projected death rate per day        

#d = (0.003 / (1 - 0.003)) * r # Projected death rate per day  

states[ , date := as.Date(date)]

# Extend dataset for each state until end of projection period
temp <- states[ , .(date = seq.Date(min(date),as.Date("2021-04-01"), by = "day")), by = state]


setkey(temp,state,date)
setkey(states,state,date)
  
states <- states[temp]

states <- Rtable[states]
states <- Ttable[states]

# Extend population data

setnames(states,"POPESTIMATE2019","population")
states[ , population := na.locf(population), by = state]
states[ , STATE := na.locf(STATE), by = state]
states[ , Area := na.locf(Area), by = state]
states[ , deathRate := na.locf(deathRate), by = state]
# Make projections

project <- function(asymptomatics,cases,deaths,vulnerablePopulation,deathRate,impliedT,factor)  {

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
    
    projectedNewAsymptomatics[i] <- f * factor * impliedT[i] * (projectedVulnerablePopulation[i-1]) * projectedAsymptomatics[i-1]
    projectedAsymptomatics[i] <- projectedAsymptomatics[i-1] * (1-r - w) + projectedNewAsymptomatics[i]
    
    projectedNewCases[i] <- (w + (1-f) * factor * impliedT[i] * (projectedVulnerablePopulation[i-1])) * projectedAsymptomatics[i-1]
    projectedCases[i] <- projectedCases[i-1] * (1-r) + projectedNewCases[i] - projectedNewDeaths[i-1]
    
    projectedNewDeaths[i] <- deathRate[i] * projectedCases[i]
    projectedDeaths[i] <- projectedDeaths[i-1] + projectedNewDeaths[i]

    projectedVulnerablePopulation[i] <- projectedVulnerablePopulation[i-1] - projectedNewAsymptomatics[i] - projectedNewCases[i]
        
  }
  
  return(list(projectedNewAsymptomatics,projectedAsymptomatics,projectedNewCases,projectedCases,projectedNewDeaths,projectedDeaths,projectedVulnerablePopulation))
  
}

states[ , c("projectedNewAsymptomatics","projectedAsymptomatics","projectedNewCases","projectedCases","projectedNewDeaths","projectedDeaths","projectedVulnerablePopulation") := project(asymptomatics,cases_ongoing,deaths,vulnerablePopulation,deathRate,impliedT_last,factor), by = state ]
states[ , c("projectedNewAsymptomatics_10","projectedAsymptomatics_10","projectedNewCases_10","projectedCases_10","projectedNewDeaths_10","projectedDeaths_10","projectedVulnerablePopulation_10") := project(asymptomatics,cases_ongoing,deaths,vulnerablePopulation,deathRate,impliedT_last,1), by = state ]
states[ , c("projectedNewAsymptomatics_07","projectedAsymptomatics_07","projectedNewCases_07","projectedCases_07","projectedNewDeaths_07","projectedDeaths_07","projectedVulnerablePopulation_07") := project(asymptomatics,cases_ongoing,deaths,vulnerablePopulation,deathRate,impliedT_last,0.7), by = state ]
states[ , c("projectedNewAsymptomatics_05","projectedAsymptomatics_05","projectedNewCases_05","projectedCases_05","projectedNewDeaths_05","projectedDeaths_05","projectedVulnerablePopulation_05") := project(asymptomatics,cases_ongoing,deaths,vulnerablePopulation,deathRate,impliedT_last,0.5), by = state ]
states[ , c("projectedNewAsymptomatics_03","projectedAsymptomatics_03","projectedNewCases_03","projectedCases_03","projectedNewDeaths_03","projectedDeaths_03","projectedVulnerablePopulation_03") := project(asymptomatics,cases_ongoing,deaths,vulnerablePopulation,deathRate,impliedT_last,0.3), by = state ]


# Make plots

states[ projectedCases == 0, projectedCases := NA]
states[ projectedNewCases == 0, projectedNewCases := NA]
states[ projectedAsymptomatics == 0, projectedAsymptomatics := NA]
states[ projectedDeaths == 0, projectedDeaths := NA]
states[ projectedVulnerablePopulation == 0, projectedVulnerablePopulation := NA]

states[ cases_ongoing == 0, cases_ongoing := NA]
states[ newCasesRaw == 0, newCasesRaw := NA]
states[ asymptomatics == 0, asymptomatics := NA]
states[ deaths == 0, deaths := NA]
states[, cumulativeInfected := population - vulnerablePopulation]

states[ , projectedCumulativeInfected := population - projectedVulnerablePopulation]
states[ , projectedCumulativeInected_10 := population - projectedVulnerablePopulation_10]
states[ , projectedCumulativeInected_07 := population - projectedVulnerablePopulation_07]
states[ , projectedCumulativeInected_05 := population - projectedVulnerablePopulation_05]
states[ , projectedCumulativeInected_03 := population - projectedVulnerablePopulation_03]

#states[ projectedDeaths == 0.0, deaths := NA]

UStotals <- states[ , .(projectedCases = sum(na.omit(projectedCases)), projectedAsymptomatics = sum(na.omit(projectedAsymptomatics)), projectedDeaths = sum(na.omit(projectedDeaths)),
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
geom_line(aes(y = log(newCasesRaw), linetype = "New cases (actual, estimated)")) + 
  geom_line(aes(y = log(projectedNewCases), linetype = "New cases (projected)")) + 
  facet_wrap(~state) + 
  theme(axis.text.x = element_text(angle = 90))

ggsave("US_states_corona_newcases_projections.pdf",plot = last_plot(), width = 18, height = 13, units = "in")

ggplot(states, aes(x = date)) + 
  geom_line(aes(y = log(cumulativeInfected), linetype = "Cumulative infected (actual, estimated)")) + 
  geom_line(aes(y = log(projectedCumulativeInfected), linetype = "Cumulative infected (projected)")) + 
  facet_wrap(~state) + 
  theme(axis.text.x = element_text(angle = 90))

ggsave("US_states_corona_infected_projections.pdf",plot = last_plot(), width = 18, height = 13, units = "in")

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
  stat_density(geom = "line", color = "black", linetype = "dashed", size = 1)

ggsave("necessaryR_adjustment.pdf", plot = last_plot(), width = 10, height = 6, units = "in")



## Make US totals with different factors

UStotals <- states[STATE <= 50 , .(projectedCases_10 = sum(na.omit(projectedCases_10)),projectedAsymptomatics_10 = sum(na.omit(projectedAsymptomatics_10)), projectedDeaths_10 = sum(na.omit(projectedDeaths_10)),
                        projectedVulnerablePopulation_10 = sum(projectedVulnerablePopulation_10),
                        cases_ongoing = sum(na.omit(cases_ongoing)), asymptomatics = sum(na.omit(asymptomatics)), deaths = sum(na.omit(deaths)), vulnerablePopulation = sum(vulnerablePopulation),
                        population = sum(na.omit(population)),
                        projectedCases_07 = sum(na.omit(projectedCases_07)),projectedAsymptomatics_07 = sum(na.omit(projectedAsymptomatics_07)), projectedDeaths_07 = sum(na.omit(projectedDeaths_07)),
                        projectedVulnerablePopulation_07 = sum(projectedVulnerablePopulation_07),
                        projectedCases_05 = sum(na.omit(projectedCases_05)),projectedAsymptomatics_05 = sum(na.omit(projectedAsymptomatics_05)), projectedDeaths_05 = sum(na.omit(projectedDeaths_05)),
                        projectedVulnerablePopulation_05 = sum(projectedVulnerablePopulation_05),
                        projectedCases_03 = sum(na.omit(projectedCases_03)),projectedAsymptomatics_03 = sum(na.omit(projectedAsymptomatics_03)), projectedDeaths_03 = sum(na.omit(projectedDeaths_03)),
                        projectedVulnerablePopulation_03 = sum(projectedVulnerablePopulation_03)), by = date]

UStotals[ cases_ongoing == 0, cases_ongoing := NA]
UStotals[ asymptomatics == 0, asymptomatics := NA]
UStotals[ deaths == 0, deaths := NA]
UStotals[ vulnerablePopulation == 0, vulnerablePopulation := NA]
UStotals[, cumulativeInfected := population - vulnerablePopulation]

for (string in c("_10","_07","_05","_03")) {
  
  UStotals[ get(paste0("projectedCases",string)) == 0, (paste0("projectedCases",string)) := NA]
  UStotals[ get(paste0("projectedAsymptomatics",string)) == 0, (paste0("projectedAsymptomatics",string)) := NA]
  UStotals[ get(paste0("projectedDeaths",string)) == 0, (paste0("projectedDeaths",string)) := NA]
  UStotals[ get(paste0("projectedVulnerablePopulation",string)) == 0, (paste0("projectedVulnerablePopulation",string)) := NA]
  
  UStotals[ , (paste0("projectedCumulativeInfected",string)) := population - get(paste0("projectedVulnerablePopulation",string))]

}

UStotals <- melt(UStotals, id.vars = "date",
                 measure.vars = c("cases_ongoing","asymptomatics","deaths","vulnerablePopulation","cumulativeInfected",
                                  "projectedCases_10","projectedAsymptomatics_10","projectedDeaths_10","projectedVulnerablePopulation_10","projectedCumulativeInfected_10",
                                  "projectedCases_07","projectedAsymptomatics_07","projectedDeaths_07","projectedVulnerablePopulation_07","projectedCumulativeInfected_07",
                                  "projectedCases_05","projectedAsymptomatics_05","projectedDeaths_05","projectedVulnerablePopulation_05","projectedCumulativeInfected_05",
                                  "projectedCases_03","projectedAsymptomatics_03","projectedDeaths_03","projectedVulnerablePopulation_03","projectedCumulativeInfected_03"))

casesVariables <- c("cases_ongoing","projectedCases_10","projectedCases_07","projectedCases_05","projectedCases_03")
asymptomaticsVariables <- c("asymptomatics","projectedAsymptomatics_10","projectedAsymptomatics_07","projectedAsymptomatics_05","projectedAsymptomatics_03")
deathsVariables <- c("deaths","projectedDeaths_10","projectedDeaths_07","projectedDeaths_05","projectedDeaths_03")
vulnerablePopulationVariables <- c("vulnerablePopulation","projectedVulnerablePopulation_10","projectedVulnerablePopulation_07","projectedVulnerablePopulation_05","projectedVulnerablePopulation_03")
cumulativeInfectedVariables <- c("cumulativeInfected","projectedCumulativeInfected_10","projectedCumulativeInfected_07","projectedCumulativeInfected_05","projectedCumulativeInfected_03")

variables_raw <- c("cases_ongoing","asymptomatics","deaths","vulnerablePopulation","cumulativeInfected")
variables_10 <- c("projectedCases_10","projectedAsymptomatics_10","projectedDeaths_10","projectedVulnerablePopulation_10","projectedCumulativeInfected_10")
variables_07 <- c("projectedCases_07","projectedAsymptomatics_07","projectedDeaths_07","projectedVulnerablePopulation_07","projectedCumulativeInfected_07")
variables_05 <- c("projectedCases_05","projectedAsymptomatics_05","projectedDeaths_05","projectedVulnerablePopulation_05","projectedCumulativeInfected_05")
variables_03 <- c("projectedCases_03","projectedAsymptomatics_03","projectedDeaths_03","projectedVulnerablePopulation_03","projectedCumulativeInfected_03")

UStotals[ variable %in% casesVariables, variableType := "Cases"]
UStotals[ variable %in% asymptomaticsVariables, variableType := "Asymptomatics"]
UStotals[ variable %in% deathsVariables, variableType := "Deaths"]
UStotals[ variable %in% vulnerablePopulationVariables, variableType := "Vulnerable Population"]
UStotals[ variable %in% cumulativeInfectedVariables, variableType := "Cumulative infected"]

UStotals[ variable %in% variables_raw, projectionType := "Data"]
UStotals[ variable %in% variables_10, projectionType := "Projection: 100% current transmission"]
UStotals[ variable %in% variables_07, projectionType := "Projection: 70% current transmission"]
UStotals[ variable %in% variables_05, projectionType := "Projection: 50% current transmission"]
UStotals[ variable %in% variables_03, projectionType := "Projection: 30% current transmission"]

p1 <- ggplot(UStotals[variable %in% casesVariables], aes(x = date, y = log(value), linetype = variable)) + 
  geom_line() + 
  labs(title = "Cases")

ggsave("US_corona_cases_projections.pdf", plot = p1, width = 10, height = 6, units = "in")

p2 <- ggplot(UStotals[variable %in% asymptomaticsVariables], aes(x = date, y = log(value), linetype = variable)) + 
  geom_line()

ggsave("US_corona_asymptomatics_projections.pdf", plot = p2, width = 10, height = 6, units = "in")

ggplot(UStotals[variable %in% vulnerablePopulationVariables], aes(x = date, y = log(value), linetype = variable)) + 
  geom_line()

p3 <- ggplot(UStotals[variable %in% cumulativeInfectedVariables], aes(x = date, y = log(value), linetype = variable)) + 
  geom_line()
  
ggsave("US_corona_infections_projections.pdf", plot = p3, width = 10, height = 6, units = "in")


p4 <- ggplot(UStotals[variable %in% deathsVariables], aes(x = date, y = log(value), linetype = variable)) + 
  geom_line()

ggsave("US_corona_deaths_projections.pdf", plot = p4, width = 10, height = 6, units = "in")

UStotals[ , projectionType := factor(projectionType, levels = c("Data","Projection: 100% current transmission","Projection: 70% current transmission",
                                                                "Projection: 50% current transmission","Projection: 30% current transmission"))]

ggplot(UStotals[ variable %in% c(casesVariables,asymptomaticsVariables,cumulativeInfectedVariables,deathsVariables)], aes(x = date, y = log(value), linetype = factor(projectionType))) + 
  geom_line() +  
  labs(title = "US total coronavirus projections",
       subtitle = "Four lockdown scenarios, logarithmic terms") + 
  xlab("Date") + 
  ylab("Log(# of people)") + 
  scale_y_continuous(trans = "exp") + 
  ylim(0,NA) + 
  #geom_point(size = 0.0001, aes(shape = projectionType)) +   
  facet_wrap(~variableType) + 
  theme(axis.text.x = element_text(angle = 90))

ggsave("US_corona_all_projections_logterms.pdf", plot = last_plot(), width = 12, height = 8, units = "in")

ggplot(UStotals[ variable %in% c(casesVariables,asymptomaticsVariables,cumulativeInfectedVariables,deathsVariables)], aes(x = date, y = value, linetype = factor(projectionType))) + 
  geom_line() +  
  labs(title = "US total coronavirus projections",
       subtitle = "Four lockdown scenarios, logarithmic terms") + 
  xlab("Date") + 
  ylab("# of people") + 
  #geom_point(size = 0.0001, aes(shape = projectionType)) + 
  facet_wrap(~variableType) + 
  theme(axis.text.x = element_text(angle = 90))

ggsave("US_corona_all_projections_raw.pdf", plot = last_plot(), width = 12, height = 8, units = "in")
