#PCA of fungal climate niches.
rm(list=ls())
library(gridExtra)
library(factoextra)

#load data.----
amf <- read.csv('data/threshold_0.9_AMF_reducedFeatures.csv')
emf <- read.csv('data/threshold_0.9_EMF_reducedFeatures.csv')
tax <- read.csv('data/Taxonomy_2020_02_01.csv')

#colors
amf.col <- 'deepskyblue4'
emf.col <- 'chartreuse4'

#subset to climate niche traits for climate pca data subsets.----
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

#plot PCA.----
#pick colors.
col.lab <- c(rep(amf.col, nrow(amf_pca.dat)), rep(emf.col, nrow(emf_pca.dat)))
#group labels.
lab <- c(rep('amf'  , nrow(amf_pca.dat)), rep('emf'  , nrow(emf_pca.dat)))

#axis labels.
p.var.1 <- round(summary(all.PCA)$importance[2,1]*100,1) #grab proportion variance explained PC1.
p.var.2 <- round(summary(all.PCA)$importance[2,2]*100,1) #grab proportion variance explained PC2.

#drop plot.
all.PCA.plot <- fviz_pca_biplot(all.PCA, geom='point', repel = T, pointshape =16, pointsize = 2,
                                habillage = as.factor(lab), palette = c(amf.col, emf.col),       #color by functional group.
                                xlab = paste0('Climate PC1 (',p.var.1,'% variance explained)'),  #x-axis label.
                                ylab = paste0('Climate PC2 (',p.var.2,'% variance explained)'),  #y-axis label.
                               title = NULL) +
  labs(tag = 'A') + 
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +    #drop gridlines
  theme(axis.text=element_text(size=10)) +    #increase axis text size.
  theme(axis.line = element_line(colour = "black")) # add axis line

all.PCA.plot
