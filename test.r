

data = Theoph

# get subset of Theoph that is equal to 1
data.test = subset(data, Subject == 1 )

# plot Time x conc
plot(data.test$Time,data.test$conc)

# necessary library for nlmeODE
library(nlmeODE)

# create groupedData object that separates the data per Subject
data.grp = groupedData( conc ~ Time | Subject, data = data )

# ploting the data.grp and see all subjects
plot(data.grp)