#Sample Code of Cusp Model Analysis in R
#setwd("~/Dropbox/Other papers/Cusp Hypothesis/Submission")
setwd("~/Dropbox (MIT)/For George")
###First load packages###
#Use a mirror and load cusp package in R through Cran
library(cusp)
library(dplyr)
library(foreign)


####Read and merge longitudinal data###
cusp = read.spss("read one or more4.sav", to.data.frame=TRUE)
imag<-read.csv("ImagingKids_FinalData_T1_T4_Nov2016.csv",na.strings=c("NA","NaN", "9999","8888","999","7777"))
screen<-read.csv("Wave1-5_AllDataMaster_Oct3_repaired.csv",na.strings=c("NA","NaN", "9999","8888","999","7777"))


long=merge(cusp, screen,all=TRUE,'READ_ID')
long=merge(x = long, y = imag[ , c("READ_ID", "W3WIss_T4")], by = "READ_ID", all.x=TRUE)
#long=merge(cusp, imag,all=TRUE,'READ_ID')


#long$WR<-(rowMeans(long[c('W3WAss_T4','W3WIss_T4','TWSss_T4','TPDEss_T4')], na.rm=TRUE))

#long<-long%>%filter(W3WIraw>0|WRWIraw>0)
#screen$WIDraw<-as.numeric(screen$WIDraw)
#long<-long%>%filter(WIDraw>0)


####Run single variable models Using Cusp Package in R###

long9<-na.omit(long%>%dplyr::select('WRWIss.y','RANLss.y','W3WIss_T4','KBITss.y'))
long10<-na.omit(long%>%dplyr::select('WIDraw','RANOss.y','W3WIss_T4','KBITss.y'))
long11<-na.omit(long%>%dplyr::select('WIDraw','RANCss.y','W3WIss_T4','KBITss.y'))
long12<-na.omit(long%>%dplyr::select('WIDraw','CTELss0s.y','W3WIss_T4','KBITss.y'))
long13<-na.omit(long%>%dplyr::select('WIDraw','CTBWss0s.y','W3WIss_T4','KBITss.y'))
long14<-na.omit(long%>%dplyr::select('WIDraw','CTNRss0s.y','W3WIss_T4','KBITss.y'))


#RANL (m9)
#long10$RANLss<-log2(long10$RANLss)

fit.long9 <- cusp(y ~ W3WIss_T4, alpha ~KBITss.y+WRWIss, beta ~ RANLss,
                  data = long9)
summary(fit.long9)
dev.off()
plot(fit.long9)

plot(fit.long9, what='bifurcation', box=TRUE, axes=FALSE)
cusp3d(fit.long9, B=5, n.surf=50, theta=150)


###RAN O (m10)
#long10$RANOss<-log2(long10$RANOss)

fit.long10 <- cusp(y ~ W3WIss_T4, alpha ~KBITss.y+ WRWIss, beta ~ RANOss,
                   data = long10)
summary(fit.long10)
dev.off()
plot(fit.long10)
plot(fit.long10, what='bifurcation', box=TRUE, axes=FALSE)
cusp3d(fit.long10, B=5, n.surf=50, theta=150)


###RAN C (m11)
long11$RANCss<-log2(long11$RANCss)

fit.long11 <- cusp(y ~ W3WIss_T4, alpha ~ KBITss.y+WRWIss, beta ~ RANCss,
                  data = long11)
summary(fit.long11)
dev.off()
plot(fit.long11)
plot(fit.long11, what='bifurcation', box=TRUE, axes=FALSE)
cusp3d(fit.long11, B=5, n.surf=50, theta=150)


#CTEL
fit.long12 <- cusp(y ~ W3WIss_T4, alpha ~ KBITss.y+WRWIss, beta ~ CTELss0s.x,
                  data = long12)
summary(fit.long12)
dev.off()
plot(fit.long12)
plot(fit.long12, what='bifurcation', box=TRUE, axes=FALSE)
cusp3d(fit.long12, B=5, n.surf=50, theta=150)

##CTBW
#long13$W3WIss_T42<-log(long13$W3WIss_T4)

fit.long13 <- cusp(y ~ W3WIss_T4, alpha ~ KBITss.y+WRWIss, beta ~ CTBWss0s.x,
                  data = long13)
summary(fit.long13)
dev.off()
plot(fit.long13)
plot(fit.long13, what='bifurcation', box=TRUE, axes=FALSE)
cusp3d(fit.long13, B=5, n.surf=50, theta=150)

##CTNR

fit.long14 <- cusp(y ~ W3WIss_T4, alpha ~ KBITss.y+WRWIss, beta ~ CTNRss0s.x,
                   data = long14)
summary(fit.long14)
dev.off()
plot(fit.long14)
plot(fit.long14, what='bifurcation', box=TRUE, axes=FALSE)
cusp3d(fit.long14, B=5, n.surf=50, theta=150)

#### Run Interaction models ####
long15<-na.omit(long%>%dplyr::select('WRWIss','RANLss','CTELss0s.x','W3WIss_T4','KBITss.y'))
long16<-na.omit(long%>%dplyr::select('WRWIss','RANLss','CTBWss0s.x','W3WIss_T4','KBITss.y'))
long17<-na.omit(long%>%dplyr::select('WRWIss','RANOss','CTELss0s.x','W3WIss_T4','KBITss.y'))
long18<-na.omit(long%>%dplyr::select('WRWIss','RANOss','CTBWss0s.x','W3WIss_T4','KBITss.y'))
long19<-na.omit(long%>%dplyr::select('WRWIss','RANLss','CTELss0s.x','RANOss',
                                     'CTBWss0s.x','W3WIss_T4','KBITss.y','CTNRss0s.x','RANCss'))

#Calculate z scores and interaction term
long15$CTELss0s<-scale(long15$CTELss0s.x)
long15$RANLss<-scale(long15$RANLss)
long15$int<-long15$CTELss0s*long15$RANLss

##CTEL and RANLL
fit.long15 <- cusp(y ~ W3WIss_T4, alpha ~ KBITss.y+WRWIss, beta ~ RANLss+CTELss0s.x+int,
                   data = long15)
summary(fit.long15)
dev.off()
#png("~/Dropbox/Other papers/Cusp Hypothesis/Submission/figures/model15.png", width = 28, height = 21, units = 'cm', res = 300)
plot(fit.long15) #significant model
cusp3d(fit.long15, B=5, n.surf=50, theta=150)#error with graph


##CTBW and RANL 
#Calculate z scores and interaction term
long16$CTBWss0s.x<-scale(long16$CTBWss0s.x)
long16$RANLss<-scale(long16$RANLss)
long16$int<-long16$CTBWss0s.x*long16$RANLss

fit.long16 <- cusp(y ~ W3WIss_T4, alpha ~ KBITss+WRWIss, beta ~ RANLss+CTBWss0s.x+int,
                   data = long16)
summary(fit.long16)
#dev.off()
#png("~/Dropbox/Other papers/Cusp Hypothesis/Submission/figures/model16.png", width = 28, height = 21, units = 'cm', res = 300)
plot(fit.long16) #sig
plot(fit.long16, what='bifurcation', box=TRUE, axes=FALSE)
cusp3d(fit.long16, B=5, n.surf=50, theta=150)#error with graph


##CTEL and RANO Int
long17$CTELss0s.x<-scale(long17$CTELss0s.x)
long17$RANOss<-scale(long17$RANOss)
long17$int<-long17$CTELss0s.x*long17$RANOss
fit.long17 <- cusp(y ~ W3WIss_T4, alpha ~ KBITss+WRWIss, beta ~ RANOss+CTELss0s.x+int,
                   data = long17)
summary(fit.long17)
#dev.off()
#png("~/Dropbox/Other papers/Cusp Hypothesis/Submission/figures/model17.png", width = 28, height = 21, units = 'cm', res = 300)
plot(fit.long17) #sig birufication
plot(fit.long17, what='bifurcation', box=TRUE, axes=FALSE)
cusp3d(fit.long17, B=5, n.surf=50, theta=150)#error with graph


##CTBW and RANO

long18$CTBWss0s.x<-scale(long18$CTBWss0s.x)
long18$RANOss<-scale(long18$RANOss)
long18$int<-long18$CTBWss0s.x*long18$RANOss

fit.long18 <- cusp(y ~ W3WIss_T4, alpha ~ KBITss+WRWIss, beta ~ RANOss+CTBWss0s.x+int,
                   data = long18)
summary(fit.long18)
#dev.off()
#png("~/Dropbox/Other papers/Cusp Hypothesis/Submission/figures/model18.png", width = 28, height = 21, units = 'cm', res = 300)
plot(fit.long18) #NA
plot(fit.long18, what='bifurcation', box=TRUE, axes=FALSE)
cusp3d(fit.long18, B=5, n.surf=50, theta=150)

fit.long19 <- cusp(y ~ W3WIss_T4, alpha ~ WRWIss, beta ~ RANOss+CTBWss0s.x+CTNRss0s.x+CTELss0s.x+RANLss+RANCss,
                   data = long19)
summary(fit.long19)
#dev.off()
#ng("~/Dropbox/Other papers/Cusp Hypothesis/Submission/figures/model19.png", width = 28, height = 21, units = 'cm', res = 300)
plot(fit.long19) #NA
plot(fit.long18, what='bifurcation', box=TRUE, axes=FALSE)
cusp3d(fit.long18, B=5, n.surf=50, theta=150)

