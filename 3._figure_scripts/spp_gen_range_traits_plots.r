#plotting species and genus range traits.
rm(list=ls())

#load data.----
amf <- read.csv('data/threshold_0.9_AMF_reducedFeatures.csv')
emf <- read.csv('data/threshold_0.9_EMF_reducedFeatures.csv')
tax <- read.csv('data/Taxonomy_2020_02_01.csv')

#global plot settings.----
par(mfrow = c(2,2))
#colors
amf.col <- 'deepskyblue4'
emf.col <- 'chartreuse4'

#means plot.----
x.emf <- emf$CHELSA_BIO_Annual_Mean_Temperature_mean
x.amf <- amf$CHELSA_BIO_Annual_Mean_Temperature_mean
y.emf <- emf$CHELSA_BIO_Annual_Precipitation_mean
y.amf <- amf$CHELSA_BIO_Annual_Precipitation_mean
#labels.
lab.x <- 'Temperature Niche Optimum'
lab.y <- 'Precipitation Niche Optimum'
#axis limits.
limy <- c(min(y.emf,y.amf, na.rm=T), max(y.emf,y.amf, na.rm=T))
limx <- c(min(x.emf,x.amf, na.rm=T), max(x.emf,x.amf, na.rm=T))

#plot.
plot(y.emf ~ x.emf, bty = 'l', pch = 16, xlim = limx, ylim = limy, xlab=lab.x, ylab=lab.y, col = emf.col)
par(new=T)
plot(y.amf ~ x.amf, bty = 'l', pch = 16, col = amf.col, xlim = limx, ylim = limy, xlab=NA, ylab=NA)

#Plot legend.
legend('topleft',legend = c('AMF','EMF'), 
       pch = 16, col = c(amf.col, emf.col), bty = 'n', cex = 1,
       x.intersp = .75, xpd = T, 
       horiz = F)


#extremes plot.----
#variables.
x.emf <- emf$CHELSA_BIO_Annual_Mean_Temperature_max
x.amf <- amf$CHELSA_BIO_Annual_Mean_Temperature_max
y.emf <- emf$CHELSA_BIO_Annual_Precipitation_min
y.amf <- amf$CHELSA_BIO_Annual_Precipitation_min
#labels.
lab.x <- 'Max Annual Temperature'
lab.y <- 'Minimum Annual Precipitation'
#axis limits.
limy <- c(min(y.emf,y.amf, na.rm=T), max(y.emf,y.amf, na.rm=T))
limx <- c(min(x.emf,x.amf, na.rm=T), max(x.emf,x.amf, na.rm=T))

#plot.
plot(y.emf ~ x.emf, bty = 'l', pch = 16, xlim = limx, ylim = limy, xlab=lab.x, ylab=lab.y, col = emf.col)
par(new=T)
plot(y.amf ~ x.amf, bty = 'l', pch = 16, col = amf.col, xlim = limx, ylim = limy, xlab=NA, ylab=NA)

#plot temp width vs. moisture widith.----
#x is temperature width.
x.emf <- emf$CHELSA_BIO_Annual_Mean_Temperature_max - emf$CHELSA_BIO_Annual_Mean_Temperature_min
x.amf <- amf$CHELSA_BIO_Annual_Mean_Temperature_max - amf$CHELSA_BIO_Annual_Mean_Temperature_min
#y is moisture width.
y.emf <- emf$CHELSA_BIO_Annual_Precipitation_max - emf$CHELSA_BIO_Annual_Precipitation_min
y.amf <- amf$CHELSA_BIO_Annual_Precipitation_max - amf$CHELSA_BIO_Annual_Precipitation_min

#labels.
lab.x <- 'Temperature Niche Width'
lab.y <- 'Precipitation Niche Width'
#axis limits.
limy <- c(min(y.emf,y.amf, na.rm=T), max(y.emf,y.amf, na.rm=T))
limx <- c(min(x.emf,x.amf, na.rm=T), max(x.emf,x.amf, na.rm=T))

#plot.
#plot(y.emf ~ x.emf, bty = 'l', pch = 16, xlim = limx, ylim = limy, xlab=lab.x, ylab=lab.y, col = emf.col)
#par(new=T)
#plot(y.amf ~ x.amf, bty = 'l', pch = 16, col = amf.col, xlim = limx, ylim = limy, xlab=NA, ylab=NA)


#plot temp width vs. temperature optimum.----
#x is temperature optimum
x.emf <- emf$CHELSA_BIO_Annual_Mean_Temperature_mean
x.amf <- amf$CHELSA_BIO_Annual_Mean_Temperature_mean
#y is temperature width.
y.emf <- emf$CHELSA_BIO_Annual_Mean_Temperature_max - emf$CHELSA_BIO_Annual_Mean_Temperature_min
y.amf <- amf$CHELSA_BIO_Annual_Mean_Temperature_max - amf$CHELSA_BIO_Annual_Mean_Temperature_min

#labels.
lab.y <- 'Temperature Niche Width'
lab.x <- 'Temperature Niche Optimum'
#axis limits.
limy <- c(min(y.emf,y.amf, na.rm=T), max(y.emf,y.amf, na.rm=T))
limx <- c(min(x.emf,x.amf, na.rm=T), max(x.emf,x.amf, na.rm=T))

#plot.
plot(y.emf ~ x.emf, bty = 'l', pch = 16, xlim = limx, ylim = limy, xlab=lab.x, ylab=lab.y, col = emf.col)
par(new=T)
plot(y.amf ~ x.amf, bty = 'l', pch = 16, col = amf.col, xlim = limx, ylim = limy, xlab=NA, ylab=NA)

#plot preiciptation width vs. precipitation optimum.----
#x is precip optimum
x.emf <- emf$CHELSA_BIO_Annual_Precipitation_mean
x.amf <- amf$CHELSA_BIO_Annual_Precipitation_mean
#y is precip width.
y.emf <- emf$CHELSA_BIO_Annual_Precipitation_max - emf$CHELSA_BIO_Annual_Precipitation_min
y.amf <- amf$CHELSA_BIO_Annual_Precipitation_max - amf$CHELSA_BIO_Annual_Precipitation_min

#labels.
lab.y <- 'Precipitation Niche Width'
lab.x <- 'Precipitation Niche Optimum'
#axis limits.
limy <- c(min(y.emf,y.amf, na.rm=T), max(y.emf,y.amf, na.rm=T))
limx <- c(min(x.emf,x.amf, na.rm=T), max(x.emf,x.amf, na.rm=T))

#plot.
plot(y.emf ~ x.emf, bty = 'l', pch = 16, xlim = limx, ylim = limy, xlab=lab.x, ylab=lab.y, col = emf.col)
par(new=T)
plot(y.amf ~ x.amf, bty = 'l', pch = 16, col = amf.col, xlim = limx, ylim = limy, xlab=NA, ylab=NA)

#plot preiciptation seasonality vs. precip minimum----
#KINDA DUMB. MORE PRECIP SEASONALITY WILL ALWAYS PUSH DOWN TOTAL PRECIP, GENERATING LOWER MIN, NO?
#I GUESS MAYBE NOT BUT WHATEVER.
#x is precip optimum
x.emf <- emf$CHELSA_BIO_Precipitation_Seasonality_mean
x.amf <- amf$CHELSA_BIO_Precipitation_Seasonality_mean
#y is precip width.
y.emf <- emf$CHELSA_BIO_Annual_Precipitation_min
y.amf <- amf$CHELSA_BIO_Annual_Precipitation_min

#labels.
lab.y <- 'min MAP'
lab.x <- 'Precipitation Seasonality'
#axis limits.
limy <- c(min(y.emf,y.amf, na.rm=T), max(y.emf,y.amf, na.rm=T))
limx <- c(min(x.emf,x.amf, na.rm=T), max(x.emf,x.amf, na.rm=T))

#plot.
#plot(y.emf ~ x.emf, bty = 'l', pch = 16, xlim = limx, ylim = limy, xlab=lab.x, ylab=lab.y, col = emf.col)
#par(new=T)
#plot(y.amf ~ x.amf, bty = 'l', pch = 16, col = amf.col, xlim = limx, ylim = limy, xlab=NA, ylab=NA)


#end plot.----
#dev.off()
