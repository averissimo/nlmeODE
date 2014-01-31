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

# setting data
data <- Indometh

# create groupedData object that separates the data per Subject
data.grp = groupedData( conc ~ time | Subject, data = data )

# get subset of Theoph that is equal to 1
data.test = subset(data, Subject == 1 )

# plot Time x conc
plot(data.test$time,data.test$conc)

# ploting the data.grp and see all subjects
labels = list(x = "Time since drug administration", y = "Theophylline serum concentration")
units  = list(x = "(hr)", y = "(mg/L)")
plot(data.grp, xlab = paste(labels$x, units$x), ylab = paste(labels$y, units$y) )

# setting 1st compartiment of diff. equation (1 of 1)
TwoComp <- list( DiffEq = list( dy1dt = ~ -(k12 + k10) * y1 +  k21 * y2,
                                dy2dt = ~ -k21 * y2 + k12 * y1 ),
                 ObsEq  = list( 
                   c1 = ~ y1 ,
                   c2 = ~ 0),
                 Parms  = c("k12", "k21", "k10", "start"),
                 States = c("y1", "y2"),
                 Init   = list("start",0))

# create nlmeODE model
model = nlmeODE( TwoComp, data.grp, 
                 LogParms = T, JAC = T, SEQ = F, rtol = 0.01, atol = 0.01)

# apply nlme to model with random effects:
#  ka + ke + CL (all variables)
data.nlme.start_k12_k10 <- nlme( conc ~ model( k12, k21, k10, start, time, Subject ),
                   data    = data.grp,
                   fixed   = k12 + k21 + k10 + start ~ 1,
                   random  = pdDiag( start + k12 + k10 ~ 1 ),
                   start   = c( k12 = 0.05, k21 = -0.5, k10 = -0.10, start = 0.7 ),
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


