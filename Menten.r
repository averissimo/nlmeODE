# 
# Non-linear mixed-effects pharmacokinetic/pharmacodynamic modelling in NLME using differential equations,
#
# Christoffer W. Tornøe, Henrik Agersø, E.Niclas Jonsson, Henrik Madsen, Henrik A. Nielsen
#  Computer Methods and Programs in Biomedicine, Volume 76, Issue 1, October 2004, Pages 31-40, ISSN 0169-2607
#   http://dx.doi.org/10.1016/j.cmpb.2004.01.001.
# 
# Adapted from paper by André Veríssimo (IDMEC/Tecnico Lisboa, PT)
# http://github.com/averissimo
#
# tested with: R version 3.0.2 (2013-09-25) -- "Frisbee Sailing"

# clear all variables
rm(list=ls())

# necessary library for nlmeODE
library(deSolve)
library(nlme)
library(nlmeODE)

# setting data (sample data taken from: http://www.erithacus.co.uk/grafit/grafit_demonstration.htm )
data <- data.frame(Time = c(1, 2, 3, 4, 5, 6), conc = c(2.8, 4.05,4.9, 5.5, 6, 6.3), Subject = rep(1,6) )

# create groupedData object that separates the data per Subject
data.grp = groupedData( conc ~ Time | Subject, data = data )

# get subset of Theoph that is equal to 1
data.test = subset(data, Subject == 1 )

# plot Time x conc
plot(data.test$Time,data.test$conc)

# ploting the data.grp and see all subjects
labels = list(x = "Substract", y = "Rate")
units  = list(x = "", y = "")
plot(data.grp, xlab = paste(labels$x, units$x), ylab = paste(labels$y, units$y) )

# setting 1st compartiment of diff. equation (1 of 1)
OneComp <- list( DiffEq = list( dy1dt = ~ Vmax * y1 / Km + y1),
                 ObsEq  = list( 
                   c1 = ~ 0),
                 Parms  = c("Vmax", "Km"),
                 States = c("y1"),
                 Init   = list(0))

# create nlmeODE model
model = nlmeODE( OneComp, data.grp, 
                 LogParms = T, JAC = T, SEQ = F, rtol = 0.01, atol = 0.01)

# apply nlme to model with random effects:
#  ka + ke + CL (all variables)
data.nlme.Vmax_Km <- nlme( conc ~ model( Vmax, Km, Time, Subject ),
                   data    = data.grp,
                   fixed   = Vmax + Km ~ 1,
                   random  = pdDiag( Vmax ~ 1 ),
                   start   = c( Vmax = 0, Km = 2.1),
                   control = list( returnObject = T, msVerbose = T ),
                   verbose = T)

# summary for the results of nlme
summary( data.nlme.start_k12_k10 )
# plot of the results, with population and subject fitting
plot(augPred(data.nlme.start_k12_k10, level=0:1))

# updates nlme using random effects:
#  ka + CL
data.nlme.start_k12 = update(data.nlme.start_k12_k10, random = pdDiag(start + k12 ~ 1))
# summary for the results of nlme
summary( data.nlme.start_k12 )
# plot of the results, with population and subject fitting
plot(augPred(data.nlme.start_k12, level=0:1))


