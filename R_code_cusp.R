#Sample Code of Cusp Model Analysis in R

###First load package###
#Use a mirror and load cusp package in R through Cran

####Then Read data###
read.table("ran2.dat", header = TRUE, row.names = 1)
ran2 = read.table("ran2.dat", header = TRUE, row.names = 1)
attach(ran2)
save(ran2, file="ran2.rda")

####Run Model Using Cusp Package in R###
data(ran2)
fit.ran2 <- cusp(y ~ WRWIss, alpha ~ WRLIss, beta ~ RANLss,
                 data = ran2)
summary(fit.ran2)
plot(fit.ran2)
plot(fit.ran2, what='bifurcation', box=TRUE, axes=FALSE)
cusp3d(fit.ran2, B=5, n.surf=50, theta=150)
