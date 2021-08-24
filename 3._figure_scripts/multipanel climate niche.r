#multi-panel Climate Niche plot.
rm(list=ls())
library(ggplot2)
library(factoextra)
library(gridExtra)
library(mgcv)
library(rlang)

#set output path.
#output.path <- figures/climate_niche.jpeg

#load data.----
amf <- read.csv('data/threshold_0.9_AMF_reducedFeatures.csv')
emf <- read.csv('data/threshold_0.9_EMF_reducedFeatures.csv')
tax <- read.csv('data/Taxonomy_2020_02_01.csv')

#colors
amf.col <- 'deepskyblue4'
emf.col <- 'chartreuse4'

#convert temperatures to degrees celcius by dividing by 10.----
amf$CHELSA_BIO_Annual_Mean_Temperature_max  <- amf$CHELSA_BIO_Annual_Mean_Temperature_max /10
amf$CHELSA_BIO_Annual_Mean_Temperature_mean <- amf$CHELSA_BIO_Annual_Mean_Temperature_mean/10
amf$CHELSA_BIO_Annual_Mean_Temperature_min  <- amf$CHELSA_BIO_Annual_Mean_Temperature_min /10
emf$CHELSA_BIO_Annual_Mean_Temperature_max  <- emf$CHELSA_BIO_Annual_Mean_Temperature_max /10
emf$CHELSA_BIO_Annual_Mean_Temperature_mean <- emf$CHELSA_BIO_Annual_Mean_Temperature_mean/10
emf$CHELSA_BIO_Annual_Mean_Temperature_min  <- emf$CHELSA_BIO_Annual_Mean_Temperature_min /10

#Setup PCA data.----
amf_pca.dat <- amf[,grep('CHELSA', colnames(amf))]
rownames(amf_pca.dat) <- amf$Taxon
amf_pca.dat <- amf_pca.dat[complete.cases(amf_pca.dat),]
emf_pca.dat <- emf[,grep('CHELSA', colnames(emf))]
rownames(emf_pca.dat) <- emf$Taxon
emf_pca.dat <- emf_pca.dat[complete.cases(emf_pca.dat),]
all_pca.dat <- rbind(amf_pca.dat, emf_pca.dat)

#calculate some new variables.
all_pca.dat$MAT_width <- all_pca.dat$CHELSA_BIO_Annual_Mean_Temperature_max - all_pca.dat$CHELSA_BIO_Annual_Mean_Temperature_min
all_pca.dat$MAP_width <- all_pca.dat$CHELSA_BIO_Annual_Precipitation_max - all_pca.dat$CHELSA_BIO_Annual_Precipitation_min  

#relabel variables.
new.lab <- c('max MAT',
             'mean MAT',
             'median MAT',
             'min MAT',
             'MAT standard deviation',
             'max MAP',
             'mean MAP',
             'median MAP',
             'min MAP',
             'MAP standard deviation',
             'max Temperature Warmest Month',
             'mean Temperature Warmest Month',
             'median Temperature Warmest Month',
             'min Temperature Warmest Month',
             'Temperature Warmest Month standard deviation',
             'max Precipitation Seasonality',
             'mean Precipitation Seasonality',
             'median Precipitation Seasonality',
             'min Precipitation Seasonality',
             'Precipitation Seasonality standard deviation',
             'MAT niche width',
             'MAP niche width'
)
colnames(all_pca.dat) <- new.lab

#subset to ones you like.
to_keep <- c(
  'mean Precipitation Seasonality',
  'max Precipitation Seasonality' ,
  'min Precipitation Seasonality' ,
  'max MAT' ,
  'max MAP' ,
  'mean MAP',
  'mean MAT',
  'min MAP' ,
  'min MAT' ,
  'MAT niche width',
  'MAP niche width'
)
all_pca.dat <- all_pca.dat[,colnames(all_pca.dat) %in% to_keep]


#Run PCA.----
all.PCA <- prcomp(all_pca.dat, center = T, scale = T)

#Panel 1. Means plot.----
x.emf <- emf$CHELSA_BIO_Annual_Mean_Temperature_mean
x.amf <- amf$CHELSA_BIO_Annual_Mean_Temperature_mean
y.emf <- emf$CHELSA_BIO_Annual_Precipitation_mean
y.amf <- amf$CHELSA_BIO_Annual_Precipitation_mean
y <- c(y.emf, y.amf)
x <- c(x.emf, x.amf)
lab <- c(rep('emf', length(x.emf)), rep('amf',length(x.amf)))
scatter.dat <- data.frame(x,y,lab)
scatter.dat <- scatter.dat[complete.cases(scatter.dat),]
col.lab <- c(rep(emf.col, nrow(scatter.dat[scatter.dat$lab == 'emf',])), 
             rep(amf.col, nrow(scatter.dat[scatter.dat$lab == 'amf',]))
                 )
#gam fits for AM and EM subsets.
fit.AM <- gam(y.amf ~ s(x.amf))
fit.EM <- gam(y.emf ~ s(x.emf))
fit.dat.AM <- data.frame(predict(fit.AM), na.omit(x.amf))
fit.dat.EM <- data.frame(predict(fit.EM), na.omit(x.emf))
colnames(fit.dat.AM) <- c('y','x')
colnames(fit.dat.EM) <- c('y','x')
rsq.amf <- round(summary(fit.AM)$r.sq, 2)
rsq.emf <- round(summary(fit.EM)$r.sq, 2)
lab.amf <- expression(R['AMF']^{"2"}==0.95)
lab.emf <- expression(R['EMF']^{"2"}==0.48)

#labels.
lab.x <- "Temperature Niche Optimum (°C)"
lab.y <- bquote('Precipitation Niche Optimum (mm yr'^-1*')')

panel.1 <- ggplot(scatter.dat, aes(x=x, y=y), size = 1) + 
  geom_point(colour = col.lab) + 
  theme_bw()  + #drop gray background.
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +     #drop gridlines
  xlab(lab.x) +                                                                       #x axis label.
  ylab(lab.y) +                                                                       #y axis label.
  theme(axis.line = element_line(colour = "black"), panel.border = element_blank()) + #add x-y axes, drop bounding box. 
  scale_y_continuous(expand = expansion(mult = c(.02, .02))) +                        #change where y-axis cuts off.
  scale_x_continuous(expand = expansion(mult = c(.02, .02)))  +                       #change where x-axis cuts off.
  theme(axis.text.x=element_text(size=rel(0.9))) +                                    #reduce x-axis text size.
  geom_smooth(data=fit.dat.AM, col = amf.col) +                                       #add AM spline fit.
  geom_smooth(data=fit.dat.EM, col = emf.col) +                                       #add EM spline fit.
  annotate(geom='text', x=-4, y=2700, label=lab.amf, parse = T) +
  annotate(geom='text', x=-4, y=2550, label=lab.emf, parse = T) 
#panel.1

#Panel 2. extremes plot.----
x.emf <- emf$CHELSA_BIO_Annual_Mean_Temperature_max
x.amf <- amf$CHELSA_BIO_Annual_Mean_Temperature_max
y.emf <- emf$CHELSA_BIO_Annual_Precipitation_min
y.amf <- amf$CHELSA_BIO_Annual_Precipitation_min
y <- c(y.emf, y.amf)
x <- c(x.emf, x.amf)
lab <- c(rep('emf', length(x.emf)), rep('amf',length(x.amf)))
scatter.dat <- data.frame(x,y,lab)
scatter.dat <- scatter.dat[complete.cases(scatter.dat),]
col.lab <- c(rep(emf.col, nrow(scatter.dat[scatter.dat$lab == 'emf',])), 
             rep(amf.col, nrow(scatter.dat[scatter.dat$lab == 'amf',]))
)
#fit models.
fit.AM <- gam(y.amf ~ s(x.amf))
fit.EM <- gam(y.emf ~ s(x.emf))

#labels.
lab.x <- "Maximum Annual Temperature (°C)"
lab.y <- bquote('Minimum Annual Precipitation (mm yr'^-1*')')

panel.2 <- ggplot(scatter.dat, aes(x=x, y=y), size = 1) + 
  geom_point(colour = col.lab) + 
  theme_bw()  + #drop gray background.
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +     #drop gridlines
  xlab(lab.x) +                                                                       #x axis label.
  ylab(lab.y) +                                                                       #y axis label.
  theme(axis.line = element_line(colour = "black"), panel.border = element_blank()) + #add x-y axes, drop bounding box. 
  scale_y_continuous(expand = expansion(mult = c(.02, .02))) +                        #change where y-axis cuts off.
  scale_x_continuous(expand = expansion(mult = c(.02, .02)))  +                       #change where x-axis cuts off.
  theme(axis.text.x=element_text(size=rel(0.9)))                                      #reduce x-axis text size.
#panel.2

#Panel 3. temp width vs. temperature optimum.----
#x is temperature optimum
x.emf <- emf$CHELSA_BIO_Annual_Mean_Temperature_mean
x.amf <- amf$CHELSA_BIO_Annual_Mean_Temperature_mean
#y is temperature width.
y.emf <- emf$CHELSA_BIO_Annual_Mean_Temperature_max - emf$CHELSA_BIO_Annual_Mean_Temperature_min
y.amf <- amf$CHELSA_BIO_Annual_Mean_Temperature_max - amf$CHELSA_BIO_Annual_Mean_Temperature_min
y <- c(y.emf, y.amf)
x <- c(x.emf, x.amf)
lab <- c(rep('emf', length(x.emf)), rep('amf',length(x.amf)))
scatter.dat <- data.frame(x,y,lab)
scatter.dat <- scatter.dat[complete.cases(scatter.dat),]
col.lab <- c(rep(emf.col, nrow(scatter.dat[scatter.dat$lab == 'emf',])), 
             rep(amf.col, nrow(scatter.dat[scatter.dat$lab == 'amf',]))
)
#gam fits for AM and EM subsets.
fit.AM <- gam(y.amf ~ s(x.amf))
fit.EM <- gam(y.emf ~ s(x.emf))
fit.dat.AM <- data.frame(predict(fit.AM), na.omit(x.amf))
fit.dat.EM <- data.frame(predict(fit.EM), na.omit(x.emf))
colnames(fit.dat.AM) <- c('y','x')
colnames(fit.dat.EM) <- c('y','x')
rsq.amf <- round(summary(fit.AM)$r.sq, 2)
rsq.emf <- round(summary(fit.EM)$r.sq, 2)
lab.amf <- expression(R['AMF']^{"2"}==0.63)
lab.emf <- expression(R['EMF']^{"2"}==0.50)

#labels.
lab.x <- "Temperature Niche Optimum (°C)"
lab.y <- "Temperature Niche Width (°C)"

panel.3 <- ggplot(scatter.dat, aes(x=x, y=y), size = 1) + 
  geom_point(colour = col.lab) + 
  theme_bw()  + #drop gray background.
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +     #drop gridlines
  xlab(lab.x) +                                                                       #x axis label.
  ylab(lab.y) +                                                                       #y axis label.
  theme(axis.line = element_line(colour = "black"), panel.border = element_blank()) + #add x-y axes, drop bounding box. 
  scale_y_continuous(expand = expansion(mult = c(.02, .02))) +                        #change where y-axis cuts off.
  scale_x_continuous(expand = expansion(mult = c(.02, .02)))  +                       #change where x-axis cuts off.
  theme(axis.text.x=element_text(size=rel(0.9))) +                                    #reduce x-axis text size.
  geom_smooth(data=fit.dat.AM, col = amf.col) +                                       #add AM spline fit.
  geom_smooth(data=fit.dat.EM, col = emf.col) +                                       #add EM spline fit.
  annotate(geom='text', x=-4, y= 5, label=lab.amf, parse = T) +
  annotate(geom='text', x=-4, y= 2, label=lab.emf, parse = T) 

#panel.3

#Panel 4. precip niche width vs. precip niche optimum.----
#x is precip optimum
x.emf <- emf$CHELSA_BIO_Annual_Precipitation_mean
x.amf <- amf$CHELSA_BIO_Annual_Precipitation_mean
#y is precip width.
y.emf <- emf$CHELSA_BIO_Annual_Precipitation_max - emf$CHELSA_BIO_Annual_Precipitation_min
y.amf <- amf$CHELSA_BIO_Annual_Precipitation_max - amf$CHELSA_BIO_Annual_Precipitation_min
y <- c(y.emf, y.amf)
x <- c(x.emf, x.amf)
lab <- c(rep('emf', length(x.emf)), rep('amf',length(x.amf)))
scatter.dat <- data.frame(x,y,lab)
scatter.dat <- scatter.dat[complete.cases(scatter.dat),]
col.lab <- c(rep(emf.col, nrow(scatter.dat[scatter.dat$lab == 'emf',])), 
             rep(amf.col, nrow(scatter.dat[scatter.dat$lab == 'amf',]))
)
#gam fits for AM and EM subsets.
fit.AM <- gam(y.amf ~ s(x.amf))
fit.EM <- gam(y.emf ~ s(x.emf))
fit.dat.AM <- data.frame(predict(fit.AM), na.omit(x.amf))
fit.dat.EM <- data.frame(predict(fit.EM), na.omit(x.emf))
colnames(fit.dat.AM) <- c('y','x')
colnames(fit.dat.EM) <- c('y','x')
rsq.amf <- round(summary(fit.AM)$r.sq, 2)
rsq.emf <- round(summary(fit.EM)$r.sq, 2)
lab.amf <- expression(R['AMF']^{"2"}==0.53)
lab.emf <- expression(R['EMF']^{"2"}==0.25)

#labels.
lab.x <- bquote('Precipitation Niche Optimum (mm yr'^-1*')')
lab.y <- bquote('Precipitation Niche Width (mm yr'^-1*')')

panel.4 <- ggplot(scatter.dat, aes(x=x, y=y), size = 1) + 
  geom_point(colour = col.lab) + 
  theme_bw()  + #drop gray background.
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +     #drop gridlines
  xlab(lab.x) +                                                                       #x axis label.
  ylab(lab.y) +                                                                       #y axis label.
  theme(axis.line = element_line(colour = "black"), panel.border = element_blank()) + #add x-y axes, drop bounding box. 
  scale_y_continuous(expand = expansion(mult = c(.02, .02))) +                        #change where y-axis cuts off.
  scale_x_continuous(expand = expansion(mult = c(.02, .02)))  +                       #change where x-axis cuts off.
  theme(axis.text.x=element_text(size=rel(0.9))) +                                    #reduce x-axis text size.
  geom_smooth(data=fit.dat.AM, col = amf.col) +                                       #add AM spline fit.
  geom_smooth(data=fit.dat.EM, col = emf.col) +                                       #add EM spline fit.
  annotate(geom='text', x=2500, y=1000, label=lab.amf, parse = T) +
  annotate(geom='text', x=2500, y= 500, label=lab.emf, parse = T) 


#Panel 5. Climate PCA plot.----
col.lab <- c(rep(amf.col, nrow(amf_pca.dat)), rep(emf.col, nrow(emf_pca.dat)))
#group labels.
lab <- c(rep('amf'  , nrow(amf_pca.dat)), rep('emf'  , nrow(emf_pca.dat)))

#axis labels.
p.var.1 <- round(summary(all.PCA)$importance[2,1]*100,1) #grab proportion variance explained PC1.
p.var.2 <- round(summary(all.PCA)$importance[2,2]*100,1) #grab proportion variance explained PC2.

#define plot.
panel.5      <- fviz_pca_biplot(all.PCA, geom='point', repel = T, pointshape =16, pointsize = 2, xlim=c(-10,10), ylim=c(-8,8),
                                habillage = as.factor(lab), palette = c(amf.col, emf.col),       #color by functional group.
                                xlab = paste0('Climate PC1 (',p.var.1,'% variance explained)'),  #x-axis label.
                                ylab = paste0('Climate PC2 (',p.var.2,'% variance explained)'),  #y-axis label.
                                title = NULL) +
  labs(tag = 'A') + 
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +    #drop gridlines
  theme(axis.text=element_text(size=10)) +    #increase axis text size.
  theme(axis.line = element_line(colour = "black")) # add axis line

#Setup plot grid and save.----
#save line.
jpeg('climate_niche_multipanel.jpeg', width=14, height= 7, units='in', res=300)
#drop panels.
grid.arrange(panel.1, panel.2, panel.3, panel.4, panel.5,
             widths = c(3,3,6),
             heights = c(3,3),
             layout_matrix = rbind(c(1,2,5),
                                   c(3,4,5)
                                  )
            )
#end plot.
dev.off()
