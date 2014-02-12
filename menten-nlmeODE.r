
# clear all variables
rm(list=ls())

# load necessary library
library(deSolve)
library(nlmeODE)

# Michaelis-Menten model
#
# parms: Vmax, Km
# States: S

# Parameters
parameters <- c( Vmax = 2.4, Km = 3.3 )

# States
state <- c( S = 8 )

# definition of function michaelis-menten, we called it menten
menten <- function(t, state, parameters) {
  with(as.list(c(state,parameters)), {
    #
    dS =  - Vmax * S / (Km + S)
    #
    list(c(dS))
  })
}

# time array
time_parms = list ( step = 0.5 , from = 0 , to = 5 , end = 15.1); time_parms$max_step = 1
times <- c(seq(time_parms$from                , time_parms$to , by = time_parms$step ), 
           seq(time_parms$to + time_parms$step, time_parms$end, by = time_parms$max_step ))

# generate data
out <- ode( y = state, times = times, func = menten, parms = parameters )

# plot the data
par( oma = c(0, 1, 3, 0))
plot(out, xlab = "Time", ylab="Concentration", type="p" )
mtext(outer=T, side = 3, "Michaelis-Menten model", cex = 1.5)

# randomize data

Time = c(); conc = c() ; Subject = c(); len = length(out[,"time"]);
for( x in 1:5 ) {
  Time    = c(Time   , out[,"time"])
  Subject = c(Subject, rep(x, len))
  conc    = c(conc   , out[,"S"] + rnorm(len) / (x*5))
}

# create data frame and grouped object
random_data     = data.frame(Time,conc,Subject)
random_data.grp = groupedData( conc ~ Time | Subject)
plot(random_data.grp)

# One Compartiment Models
OneComp <- list( DiffEq = list( dSdt = ~ - (Vmax + S) / (Km + S) ),
                 ObsEq  = list( c1 = ~ S ),
                 Parms  = c("Vmax","Km"),
                 States = c("S"),
                 Init   = list(8))

menten.model <- nlmeODE( OneComp, random_data.grp,
                         LogParms = F, JAC = T, SEQ = F, rtol = 0.01, atol = 0.01)

menten.model( 1, 3, c(0,1,2,3,4,5), c(1,1,1,1,1,1))

data.nlme.ka_ke_CL <- nlme( conc ~ menten.model( Vmax, Km, Time, Subject ),
                            data    = random_data.grp,
                            fixed   = Vmax + Km ~ 1,
                            random  = Vmax + Km ~ 1 | Subject,
                            start   = c( Vmax = 2.4, Km = 3.3 ),
                            control = list( returnObject = T, msVerbose = T ),
                            verbose = T)

summary(data.nlme.ka_ke_CL)

plot(augPred(data.nlme.ka_ke_CL, level=0:1))
