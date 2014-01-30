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
data <- Theoph
# Why are these needed??!?
data$Dose[data$Time!=0] <- 0    # necessary to speed nlme and also to be able to plot(augPred(...))
data$Cmt <- rep(1,dim(data)[1]) # without this optional data column, nlme fails!

# create groupedData object that separates the data per Subject
data.grp = groupedData( conc ~ Time | Subject, data = data )

# get subset of Theoph that is equal to 1
data.test = subset(data, Subject == 1 )

# plot Time x conc
plot(data.test$Time,data.test$conc)

# ploting the data.grp and see all subjects
labels = list(x = "Time since drug administration", y = "Theophylline serum concentration")
units  = list(x = "(hr)", y = "(mg/L)")
plot(data.grp, xlab = paste(labels$x, units$x), ylab = paste(labels$y, units$y) )

# setting 1st compartiment of diff. equation (1 of 1)
OneComp <- list( DiffEq = list( dy1dt = ~ -ka * y1,
                                dy2dt = ~ ka * y1 - ke * y2 ),
                 ObsEq  = list(c1 = ~ 0, c2 = ~ y2/CL*ke),
                 #ObsEq  = list( SC ~ 0, Cp ~ y2 / CL * ke),
                 Parms  = c("ka", "ke", "CL"),
                 States = c("y1", "y2"),
                 Init   = list(0,0)
                 )

# create nlmeODE model
model = nlmeODE( OneComp, data, 
                 LogParms = T, JAC = T, SEQ = F, rtol = 0.01, atol = 0.01)

# apply nlme to model with random effects:
#  ka + ke + CL (all variables)
data.nlme.ka_ke_CL <- nlme( conc ~ model( ka, ke, CL, Time, Subject ),
                   data    = data,
                   fixed   = ka + ke + CL ~ 1,
                   random  = pdDiag( ka + ke + CL ~ 1 ),
                   start   = c( ka = log(1.65), ke = log(0.08), CL = log(0.05) ),
                   control = list( returnObject = T, msVerbose = T ),
                   verbose = T)

# summary for the results of nlme
summary( data.nlme.ka_ke_CL )
# plot of the results, with population and subject fitting
plot(augPred(data.nlme.ka_ke_CL, level=0:1))

# updates nlme using random effects:
#  ka + CL
data.nlme.ka_CL = update(data.nlme.ka_ke_CL, random = pdDiag(ka + CL ~ 1))
# summary for the results of nlme
summary( data.nlme.ka_CL )
# plot of the results, with population and subject fitting
plot(augPred(data.nlme.ka_CL, level=0:1))

# updates nlme using a more accurate error correction for this data
#  it is observerd that data.nlme.ka_CL residuals increase over time
#  as such weigths will reflect this
data.nlme.ka_CL.error2 <- update(data.nlme.ka_CL, weights = varConstPower( const = 0.7, power = 0.3))
# summary for the results of nlme
summary( data.nlme.ka_CL.error2 )
# plot of the results, with population and subject fitting
plot(augPred(data.nlme.ka_CL, level=0:1))

