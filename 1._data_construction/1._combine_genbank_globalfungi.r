#workup Stephanie's genbank data into something useable by Johan.
rm(list=ls())
library(data.table)

#Set output path.----
output.path <- 'data/EMF_presences_to_map.csv'

#load data.----
#genbank data.
d <- read.csv('data/EMFungiGenBank_2020_11_19.csv')
d$RDP_ID <- as.character(d$RDP_ID)

#globalfungi data.
#Unfortunately, these data are not publicly available at present. We will share with users downstream mapping ready data that can be used to replicate our analyses.
gf <- fread('/Users/colinaverill/Documents/git_projects/global_fungi_ECM/data/global_fungi_otu.csv', header=T)
gf.env <- read.csv('/Users/colinaverill/Documents/git_projects/global_fungi_ECM/data/global_fungi_env.csv')
gf.env$id.int <- as.character(gf.env$id.int)
gf <- as.data.frame(gf)
rownames(gf) <- gf$SHs
gf$SHs <- NULL


#grab sh data from RDP tag from genbank.----
sh.data <- ifelse(grepl('SH',d$RDP_ID), d$RDP_ID, NA)
sh.data <- sub(".*SH", "", sh.data)                 # Extract characters after pattern

#add in SH characters again.
out <- c()
for(i in 1:length(sh.data)){
  if( is.na(sh.data[i])){out[i] <- NA}
  if(!is.na(sh.data[i])){out[i] <- paste0('SH',sh.data[i])}
}

#pop into main dataframe.
d$UNITE_sh <- out

#grab genus_sp taxonomy from RDP tag.----
d$taxonomy <- sub("\\|.*", "", d$RDP_ID)


#drop non EM fungi.----
#only ~37k of 65k rows have an EM guild classification.
d$Guild <- as.character(d$Guild)
d <- d[d$Guild == 'EM',]


#turn global fungi into presence/absence table.----
gf[,] <- ifelse(gf > 1, 1, 0)
#How many global fungi SH's have been observed in >100 locations?
sub.100 <- rowSums(gf)
length(sub.100[sub.100 > 100]) #193 SHs observed > 1000 times.
sub.100.lab <- names(sub.100[sub.100 > 100])


#Check who has enough observations from genbank.----
tax.test <- table(d$taxonomy)
length(tax.test[tax.test>100])
sum(tax.test[tax.test>100])

 sh.test <- table(d$UNITE_sh)
length(sh.test[ sh.test>100])
sum( sh.test[ sh.test>100])
gb.sh.100 <- names(sh.test[sh.test > 100])

#presence table from global fungi.----
gf.presence <- list()
for(i in 1:length(sub.100.lab)){
  pres.abs <- data.frame(t(gf[rownames(gf) == sub.100.lab[i],]))
  pres.abs$sample.ID <- rownames(pres.abs)
  colnames(pres.abs)[1] <- 'present'
  pres.abs <- pres.abs[pres.abs$present == 1,]
  pres.abs <- merge(pres.abs, gf.env[,c('id.int','latitude.float','longitude.float')], by.x = 'sample.ID', by.y = 'id.int')
  pres.abs$present <- NULL
  pres.abs$unite_sh <- sub.100.lab[i]
  gf.presence[[i]] <- pres.abs
}
gf.presence <- do.call(rbind, gf.presence)
gf.presence$source <- 'global_fungi'
colnames(gf.presence) <- c('sample.ID','latitude','longitude','unite_sh','source')

#Make genbank presence table.----
gb.presence <- d[!is.na(d$UNITE_sh),]
gb.presence <- d[d$UNITE_sh %in% gb.sh.100,]
gb.presence <- gb.presence[,c('Latitude','Longitude','UNITE_sh')]
gb.presence$source <- 'genbank'
gb.presence$sample.ID <- NA
gb.presence <- gb.presence[,c('sample.ID','Latitude','Longitude','UNITE_sh','source')]
colnames(gb.presence) <- c('sample.ID','latitude','longitude','unite_sh','source')
gb.presence$sample.ID <- paste0('gb_',c(1:nrow(gb.presence)))
gb.presence$latitude  <- as.numeric(as.character(gb.presence$latitude ))
gb.presence$longitude <- as.numeric(as.character(gb.presence$longitude))
test <- gb.presence[complete.cases(gb.presence),]

#merge global fungi and genbank presence tables.----
output <- rbind(gf.presence, gb.presence)

#save output.----
write.csv(output, output.path)
