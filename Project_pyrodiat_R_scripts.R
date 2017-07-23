##################################################################################################
####SCHEFFERVILLE NGS: Sept 2015 work                                                            #
##################################################################################################

##In this script: Preliminary analyses of cleaned diatom data from
#Sept 2015 PANAM output. 
#Previous script: Projet_pyroeuk_R_scripts.R  
#R version: 3.1.2 (Pumpkin Helmet)

##Last update: Sept. 24, 2015
##Associated workspace: 
##Associated markdown: 
##Associated .txt of R script: 
##Github: 

##################################################################################################
##Working with diatom data to investigate patterns found with eukaryote data.
#Specifically to investigate shifts in diatom assemblages around ~1955. 

##################################################################################################

setwd("C:/Users/Winegardner/Documents/MCGILL/PhD chapters and projects/Schefferville NGS/For Sept 2015/Cleaning_data_Sept2015")

####Packages
library(vegan) #multivariate/ordination 
library(MASS) #multivariate/ NMDS
library(fossil) #paleo analysis tool
library(digest)
library(stats)
library(gclus) #clustering 
library(ggplot2) #graphs/visuals
library(ggdendro) #using ggplot to plot cluster analyses
library(gridExtra) #graphs/visuals 
library(plyr) #data manipulation/summary functions
library(reshape) #data manipulation
library(reshape2) #data manipulation 
library(rpart) #univariate regression trees
library(lme4) #linear mixed effect models 
library(AICcmodavg) #AIC
library(arm) #AIC/regression modelling 
library(mvpart) #multivariate regression trees (archived package)
library(MVPARTwrap) #MRT plotting (archived package)
library(rpart.plot) #Regression trees
library(partykit) #Regression trees
library(tree) #Regression trees

##START WITH DIAT_NOSING_BAC

#################RAW DATA WITH NO NORMALIZATION##################
##Setting up raw data for use for subsequent analyses without rarefying. 
#Going to use this data as is for now because of low minimum read numbers. 

raw<- read.table("diat_nosing_bac.csv",h=T, sep=",") ##diatom data, no singletons, resolved to genus
#Remove Sed 10 (column 11) since removed from eukaryote data. 
raw<- raw[,-11]


#################NORMY FUNCTION##################
##Function to "normalize" data by the minimum number of reads (rarefied read #).  
#Remove sample 10 (the last DAR sample as affecting the normalization, too low)

a <- read.table("diat_nosing_bac.csv",h=T, sep=",") #diatom data, no singletons, resolved to genus

#Remove column 11 (Sed10) - because removed from eukaryote data
a<- a[,-11]

#Look at colSums --> really do have an issue here as # of reads is so low. 
#Should also remove DAR5 (Sed 5)- if leave in the rest of the samples, would be normalizing
#to minimum of 55 reads ******STILL NOT PREFERABLE. 
a<- a[,-6] #so for now remove column 6 (Sed 5) 

data.frame(colnames(a))
min(colSums(a[,2:15])) #for minimum number of reads, now 55. 

#CONSIDER USING DATA RAW/ NON-RAREFIED.*****

normy = function(data, number_of_samples,number_of_reads){
  
  library(vegan)
  
  z <- list() 
  for(i in 1:number_of_samples){
    z[[i]] <- rrarefy(data[,(i+1)], sample = number_of_reads)}
  z2 <- matrix(unlist(z), ncol=number_of_samples)
  z3 <- cbind(data[,1], z2,data[, (number_of_samples+2):length(data)])
  abond <- rowSums(z3[,2:(number_of_samples+1)]) 
  z4 <- cbind(z3, abond)
  z5 <- subset(z4, z4[,length(z4)] !="0")
  z5 <- z5[,-length(z5)]
  colnames(z5) <- colnames(data)
  z5
}

b <- normy(data=a,number_of_samples=14,number_of_reads=55)
write.table(b,file="diat_nosing_14s55r.txt",sep="\t",row.names=F) #labelled with 14 samples and 55 reads.


#################RAREFACTION CURVES##################
##Create rarefaction curves to see aysmptotes for rarefied read counts. 

a<-read.table("diat_nosing_14s55r.txt",h=T,row.names=1) #uses output from normy
a1<-as.data.frame(t(a[,1:14]))
S<-specnumber(a1)
raremax<-min(rowSums(a1))
Srare<-rarefy(a1,raremax)
plot(S,Srare,xlab="Observed No.ofSpecies",ylab="Rarefied No.ofSpecies")
abline(0,1)
rarecurve(a1,step=100,sample=592,col="blue",cew=0.6)


#################METRY FUNCTION:CALCULATING RICHNESS AND DIVERSITY METRICS##################
##Function to calculate variety of diversity metrics on read data. 

metry = function(data, number_of_samples){
  library(fossil)
  library(vegan)
  
  a1 <- list()  
  for(i in 2:(number_of_samples+1))
  {a1[[i]] <- as.numeric(data[,i] > 0)}
  a2 <- matrix(unlist(a1), ncol = number_of_samples)
  a_bin <- cbind(data[,1], a2, data[,(number_of_samples+2):length(data)])
  colnames(a_bin) <- colnames(data)
  
  nb_otus <- colSums(a_bin[,2:(number_of_samples+1)])
  shannon <- diversity(t(data[,2:(number_of_samples+1)]), index="shannon")
  simpson <- diversity(t(data[,2:(number_of_samples+1)]), index="simpson")
  invsimpson <- diversity(t(data[,2:(number_of_samples+1)]), index="invsimpson")
  pielou <- shannon/log(specnumber(t(data[,2:(number_of_samples+1)])))
  
  z <- list()      
  for( i in 2:(number_of_samples+1)) { 
    z[[i]] <- chao1(data[,i])}      
  chao1 <- unlist(z)
  
  z <- list()      
  for( i in 2:(number_of_samples+1)) { 
    z[[i]] <- chao2(data[,i])}    	
  chao2 <- unlist(z)
  
  z <- list()      
  for( i in 2:(number_of_samples+1)) { 
    z[[i]] <- ACE(data[,i])}    	
  ACE <- unlist(z)
  
  z <- list()      
  for( i in 2:(number_of_samples+1)) { 
    z[[i]] <- ICE(data[,i], taxa.row=T)}      
  ICE <- unlist(z)
  
  z <- list() 
  for(i in 2:(number_of_samples+1)){
    z[[i]] <- nrow(subset(data, data[,i] =="1"))  }
  z2 <- matrix(unlist(z), ncol=number_of_samples)
  A <- colSums(a[,2:(number_of_samples+1)])
  Good <- data.frame(1 - z2/A)
  
  b <- as.data.frame(cbind(colnames(data[,2:(number_of_samples+1)]),nb_otus, shannon, simpson, invsimpson, pielou, chao1, chao2, ACE, ICE,Good))
  colnames(b) <- c("Samples", "Nb_otus", "Shannon", "Simpson", "Invsimpson", "Pielou", "Chao1", "Chao2", "ACE", "ICE", "Good's coverage")
  b
}

##Raw data 
#use raw data
raw #Sample columns are 2:16 (15 samples) 

raw.b<-metry(data=raw, number_of_samples=15)
write.table(raw.b[,1:11], file="diat_nosing_15sRAWr_metrics.txt", sep="\t", row.names=F)
#only need 1st 11 columns. 

##Rarefied data
a <- read.table("diat_nosing_14s55r.txt",h=T) #use rarefied data (output from normy) 
data.frame(colnames(a))

b<-metry(data=a,number_of_samples=14)
write.table(b[,1:11],file="diat_nosing_14s55r_metrics.txt",sep="\t",row.names=F)


#################TAXY:TAXO GROUPING##################
##Summarize OTUs by NN_Genus

##Raw data
raw #use raw from above 
data.frame(colnames(raw))

raw.long<- melt(raw) #OTUs, NN_Family and NN_Genus get used as id variables- good. 
colnames(raw.long)[4]<- 'Sample_ID'
colnames(raw.long)[5]<- 'No.Reads'

raw.cast<- dcast(raw.long, Sample_ID~NN_Genus, value.var = "No.Reads", fun.aggregate=sum)
#works- just have that weird extra 'X' row at bottom, remove. 
raw.cast<- raw.cast[-16,]
#numbers in columns represent # of reads. 

raw.transpose<- as.data.frame(t(raw.cast))
write.csv(raw.transpose, "Genusbysample_15sRAWr.csv")


##Rarefied data 
a<-read.table("diat_nosing_14s55r.txt",h=T) #rarefied data (output from normy)
data.frame(colnames(a))

a.long<- melt(a) #OTUs, NN_Family and NN_Genus get used as id variables- good. 
colnames(a.long)[4]<- 'Sample_ID'
colnames(a.long)[5]<- 'No.Reads'

a.cast<- dcast(a.long, Sample_ID~NN_Genus, value.var = "No.Reads", fun.aggregate=sum)
#works- just have that weird extra 'X' row at bottom, remove. 
a.cast<- a.cast[-15,]
#numbers in columns represent # of reads. 

a.transpose<- as.data.frame(t(a.cast))
write.csv(a.transpose, "Genusbysample_14s55r.csv")


#################HISTOGRAMS OF NN_Genus##################
#Create some summary plots in ggplot2. 

##Raw data
#Make raw.cast long again. 
rawcast.long<- melt(raw.cast, id.vars=c("Sample_ID"))
colnames(rawcast.long)[2]<- 'NN_Genus'
colnames(rawcast.long)[3]<- 'No_reads'

d<- ggplot(rawcast.long, aes(x=Sample_ID, y=No_reads, fill=NN_Genus)) + geom_bar(stat = 'identity')
d<- d + scale_fill_manual(values = c("Amphora" = "darkgoldenrod1", "Bacillaria" = "darkblue", "Bellerochea" = "pink", "Caloneis" = "khaki",
                                     "Cymbopleura" = "seagreen3", "Cylindrotheca" = "white", "Discostella" = "orangered3", "Encyonema" = "yellow2",
                                     "Fragilaria" = "turquoise4", "Haslea" = "slategrey", "Lauderia" = "olivedrab3",
                                     "Lemnicola" = "limegreen", "Melosira" = "darksalmon", "Nitzschia" = "chocolate2",
                                     "Phaeodactylum" = "cadetblue4", "Pinnularia" = "blue", "Punctastriata" = "darkgreen", "Sellaphora" = "darkred",
                                     "Staurosira" = "firebrick1", "Stenopterobia" = "blue", "Surirella" = "gray11", "Thalassiosira" = "deeppink",
                                     "Unclassified" = "black"))
d<- d + labs(x= "Sample", y= "Number of reads") + theme_bw()
d<- d + theme(axis.text.x = element_text(colour="black", angle=45, size=16))
d<- d + theme(axis.text.y = element_text(colour="black", size=16))
d<- d + theme(axis.title.x = element_text(size = rel(2), angle=00))
d<- d + theme(axis.title.y = element_text(size = rel(2), angle=90))
#Put DAR and KB on different plots so KB reads don't skew scale. 

#DAR only
rawcast.longDAR<- as.data.frame(subset(rawcast.long, Sample_ID == "Sed1_DAR1" | Sample_ID == "Sed2_DAR2" | Sample_ID == "Sed3_DAR3" | Sample_ID == "Sed4_DAR4" | Sample_ID == "Sed5_DAR5" | Sample_ID == "Sed6_DAR6" | Sample_ID == "Sed7_DAR7" | Sample_ID == "Sed8_DAR8" | Sample_ID == "Sed9_DAR9", drop=T))

d<- ggplot(rawcast.longDAR, aes(x=Sample_ID, y=No_reads, fill=NN_Genus)) + geom_bar(stat = 'identity')
d<- d + scale_fill_manual(values = c("Amphora" = "darkgoldenrod1", "Bacillaria" = "darkblue", "Bellerochea" = "pink", "Caloneis" = "khaki",
                                     "Cymbopleura" = "seagreen3", "Cylindrotheca" = "white", "Discostella" = "orangered3", "Encyonema" = "yellow2",
                                     "Fragilaria" = "turquoise4", "Haslea" = "slategrey", "Lauderia" = "olivedrab3",
                                     "Lemnicola" = "limegreen", "Melosira" = "darksalmon", "Nitzschia" = "chocolate2",
                                     "Phaeodactylum" = "cadetblue4", "Pinnularia" = "blue", "Punctastriata" = "darkgreen", "Sellaphora" = "darkred",
                                     "Staurosira" = "firebrick1", "Stenopterobia" = "blue", "Surirella" = "gray11", "Thalassiosira" = "deeppink",
                                     "Unclassified" = "black"))
d<- d + labs(x= "Sample", y= "Number of reads") + theme_bw()
d<- d + theme(axis.text.x = element_text(colour="black", angle=45, size=16))
d<- d + theme(axis.text.y = element_text(colour="black", size=16))
d<- d + theme(axis.title.x = element_text(size = rel(2), angle=00))
d<- d + theme(axis.title.y = element_text(size = rel(2), angle=90))

#KB only 
rawcast.longKB<- as.data.frame(subset(rawcast.long, Sample_ID == "Sed11_KB1" | Sample_ID == "Sed12_KB2" | Sample_ID == "Sed13_KB3" | Sample_ID == "Sed14_KB4" | Sample_ID == "Sed15_KB5" | Sample_ID == "Sed16_KB6", drop=T))

d<- ggplot(rawcast.longKB, aes(x=Sample_ID, y=No_reads, fill=NN_Genus)) + geom_bar(stat = 'identity')
d<- d + scale_fill_manual(values = c("Amphora" = "darkgoldenrod1", "Bacillaria" = "darkblue", "Bellerochea" = "pink", "Caloneis" = "khaki",
                                     "Cymbopleura" = "seagreen3", "Cylindrotheca" = "white", "Discostella" = "orangered3", "Encyonema" = "yellow2",
                                     "Fragilaria" = "turquoise4", "Haslea" = "slategrey", "Lauderia" = "olivedrab3",
                                     "Lemnicola" = "limegreen", "Melosira" = "darksalmon", "Nitzschia" = "chocolate2",
                                     "Phaeodactylum" = "cadetblue4", "Pinnularia" = "blue", "Punctastriata" = "darkgreen", "Sellaphora" = "darkred",
                                     "Staurosira" = "firebrick1", "Stenopterobia" = "blue", "Surirella" = "gray11", "Thalassiosira" = "deeppink",
                                     "Unclassified" = "black"))
d<- d + labs(x= "Sample", y= "Number of reads") + theme_bw()
d<- d + theme(axis.text.x = element_text(colour="black", angle=45, size=16))
d<- d + theme(axis.text.y = element_text(colour="black", size=16))
d<- d + theme(axis.title.x = element_text(size = rel(2), angle=00))
d<- d + theme(axis.title.y = element_text(size = rel(2), angle=90))


##Rarefied data
#Make a.cast long again. 
acast.long<- melt(a.cast, id.vars=c("Sample_ID"))
colnames(acast.long)[2]<- 'NN_Genus'
colnames(acast.long)[3]<- 'No_reads'

d<- ggplot(acast.long, aes(x=Sample_ID, y=No_reads, fill=NN_Genus)) + geom_bar(stat = 'identity')
d<- d + scale_fill_manual(values = c("Amphora" = "darkgoldenrod1", "Bacillaria" = "darkblue", "Caloneis" = "khaki",
                                     "Cymbopleura" = "seagreen3", "Discostella" = "orangered3", "Encyonema" = "yellow2",
                                     "Fragilaria" = "turquoise4", "Haslea" = "slategrey", "Lauderia" = "olivedrab3",
                                     "Lemnicola" = "limegreen", "Melosira" = "darksalmon", "Nitzschia" = "chocolate2",
                                     "Phaeodactylum" = "cadetblue4", "Pinnularia" = "blue", "Sellaphora" = "darkred",
                                     "Staurosira" = "firebrick1", "Stenopterobia" = "blue", "Surirella" = "gray11", "Thalassiosira" = "deeppink",
                                     "Unclassified" = "black"))
d<- d + labs(x= "Sample", y= "Number of reads") + theme_bw()
d<- d + theme(axis.text.x = element_text(colour="black", angle=45, size=16))
d<- d + theme(axis.text.y = element_text(colour="black", size=16))
d<- d + theme(axis.title.x = element_text(size = rel(2), angle=00))
d<- d + theme(axis.title.y = element_text(size = rel(2), angle=90))

#Try with a facet to see if easier to visualize.
d<- ggplot(acast.long, aes(x=Sample_ID, y=No_reads)) + geom_bar(stat = 'identity', fill="darkblue") + facet_wrap(~NN_Genus)
d<- d + theme(axis.text.x = element_text(colour="black", angle=90, size=10))
#d<- d + scale_fill_manual(values = c("Sed1_DAR1" = "darkblue", "Sed2_DAR2" = "darkblue", "Sed3_DAR3" = "darkblue",
 #                                    "Sed4_DAR4" = "darkblue", "Sed6_DAR7" = "darkblue", "Sed8_DAR8" = "darkblue",
 #                                   "Sed9_DAR9" = "darkblue", "Sed11_KB1" = "red", "Sed12_KB2" = "red", "Sed13_KB3" = "red",
 #                                 "Sed14_KB4" = "red", "Sed15_KB5" = "red", "Sed16_KB6" = "red"))


#################VISUALIZE DIVERSITY METRICS##################

##Raw data
div.raw<- read.csv("diat_nosing_15sRAWr_metrics_expanded.csv")

#Subset into lac Dauriat and Knob
div.raw.DAR<- as.data.frame(subset(div.raw, Lake == "DAR", drop=T))
div.raw.KB<- as.data.frame(subset(div.raw, Lake == "KB", drop=T))

#DAR plots
#Substitute y = and ylab =
#Nb_otus
#Shannon
#Simpson
#Invsimpson
#Pielou
#Chao1
#Chao2
#ACE
#ICE
#Goods_coverage
dar.raw.plot<- ggplot(div.raw.DAR, aes(x=Approx_year2, y=Nb_otus)) + geom_point(size=4, colour="darkred") + geom_path(size=0.75, colour="darkred")
dar.raw.plot<- dar.raw.plot + labs(x= "Estimated year", y= "# OTUs") + theme_bw()
dar.raw.plot<- dar.raw.plot + theme(axis.text.x = element_text(colour="black", size=16))
dar.raw.plot<- dar.raw.plot + theme(axis.text.y = element_text(colour="black", size=16))
dar.raw.plot<- dar.raw.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
dar.raw.plot<- dar.raw.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))

#KB plots 
#Substitute y = and ylab =
#Nb_otus
#Shannon
#Simpson
#Invsimpson
#Pielou
#Chao1
#Chao2
#ACE
#ICE
#Goods_coverage
kb.raw.plot<- ggplot(div.raw.KB, aes(x=Approx_year2, y=Nb_otus)) + geom_point(size=4, colour="darkred") + geom_path(size=0.75, colour="darkred")
kb.raw.plot<- kb.raw.plot + labs(x= "Estimated year", y= "# OTUs") + theme_bw()
kb.raw.plot<- kb.raw.plot + theme(axis.text.x = element_text(colour="black", size=16))
kb.raw.plot<- kb.raw.plot + theme(axis.text.y = element_text(colour="black", size=16))
kb.raw.plot<- kb.raw.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
kb.raw.plot<- kb.raw.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))


##Rarefied data 
div <- read.csv("diat_nosing_14s55r_metrics_expanded.csv")

#Subset into lac Dauriat and Knob
div.DAR<- as.data.frame(subset(div, Lake == "DAR", drop=T))
div.KB<- as.data.frame(subset(div, Lake == "KB", drop=T))

#DAR plots
#Substitute y = and ylab =
#Nb_otus
#Shannon
#Simpson
#Invsimpson
#Pielou
#Chao1
#Chao2
#ACE
#ICE
#Goods_coverage
dar.plot<- ggplot(div.DAR, aes(x=Approx_year2, y=Nb_otus)) + geom_point(size=4, colour="darkred") + geom_path(size=0.75, colour="darkred")
dar.plot<- dar.plot + labs(x= "Estimated year", y= "# OTUs") + theme_bw()
dar.plot<- dar.plot + theme(axis.text.x = element_text(colour="black", size=16))
dar.plot<- dar.plot + theme(axis.text.y = element_text(colour="black", size=16))
dar.plot<- dar.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
dar.plot<- dar.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))

#KB plots 
#Substitute y = and ylab =
#Nb_otus
#Shannon
#Simpson
#Invsimpson
#Pielou
#Chao1
#Chao2
#ACE
#ICE
#Goods_coverage
kb.plot<- ggplot(div.KB, aes(x=Approx_year2, y=Nb_otus)) + geom_point(size=4, colour="darkred") + geom_path(size=0.75, colour="darkred")
kb.plot<- kb.plot + labs(x= "Estimated year", y= "# OTUs") + theme_bw()
kb.plot<- kb.plot + theme(axis.text.x = element_text(colour="black", size=16))
kb.plot<- kb.plot + theme(axis.text.y = element_text(colour="black", size=16))
kb.plot<- kb.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
kb.plot<- kb.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))

#Note: depending on the diversity metric, will be the same across both raw and rarefied data.


#################ORDINATION ANALYSES##################
##Generate some basic ordinations in order to explore variation in
#diatom community assemblages and structure.

##Matrices
#Raw OTUs by samples
raw<- read.table("diat_nosing_bac.csv",h=T, sep=",") 
raw<- raw[,-11]#raw data (removing DAR10)
raw<- raw[,1:16] #only need OTUs and samples
raw<- melt(raw)
raw<- dcast(raw, variable~OTUs)
raw<- as.data.frame(cbind(div.raw[,2:5], raw[,2:296])) #bind with extra info (295 OTUs)
rownames(raw)<- as.character(raw[,2])

#Rarefied OTUs by samples 
a<-read.table("diat_nosing_14s55r.txt",h=T) #rarefied data
a<- a[,1:15] #only need OTUs and samples
a<- melt(a) #make long-format
a<- dcast(a, variable ~ OTUs) #cast back into transposed matrix. 
a<- as.data.frame(cbind(div[,2:5], a[,2:142])) #bind with extra info (141 OTUs)
rownames(a)<- as.character(a[,2])

#Raw diversity metrics with additional descriptive info
div.raw<- read.csv("diat_nosing_15sRAWr_metrics_expanded.csv") #diversity metrics on raw data

#Rarefied diversity metrics with additional descriptive info 
div <- read.csv("diat_nosing_14s55r_metrics_expanded.csv") #diversity metrics on rarefied data


##NMDS - sample by OTU (No. reads)
#Raw (more OTUs)
raw.nmds<-metaMDS(raw[,5:299],distance="bray") 
raw1 <- raw.nmds$points[,1]
raw2 <- raw.nmds$points[,2]
rawc<- as.data.frame(cbind(raw1, raw2))

raw.nmds.plot<- ggplot()
raw.nmds.plot<- raw.nmds.plot + geom_vline(x=0,colour="grey50")
raw.nmds.plot<- raw.nmds.plot + geom_hline(y=0,colour="grey50")
raw.nmds.plot<- raw.nmds.plot + geom_text(data = rawc, aes(x=raw1, y=raw2, label = raw$Approx_year2, colour = raw$Lake), size = 7)
raw.nmds.plot<- raw.nmds.plot + scale_colour_manual(values = c("DAR" = "darkblue", "KB" = "red2"), "Lake")
raw.nmds.plot<- raw.nmds.plot + theme_bw()
raw.nmds.plot<- raw.nmds.plot + labs(x= "NMDS 1", y= "NMDS 2")
raw.nmds.plot<- raw.nmds.plot + theme(axis.text.x = element_text(colour="black", size=16))
raw.nmds.plot<- raw.nmds.plot + theme(axis.text.y = element_text(colour="black", size=16))
raw.nmds.plot<- raw.nmds.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
raw.nmds.plot<- raw.nmds.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
raw.nmds.plot<- raw.nmds.plot + ggtitle("NMDS - OTUs (Reads (raw/non-rarefied) non-transformed)")
#Add if want vectors
raw.nmds.plot<- raw.nmds.plot + geom_path(data = rawc, aes(x = raw1, y= raw2, colour= raw$Lake), size =0.4, arrow=arrow(length= unit(0.4,"cm"), ends="first", type="open" ))

#Rarefied (less OTUs)
a.nmds<-metaMDS(a[,5:145],distance="bray") 
a1 <- a.nmds$points[,1]
a2 <- a.nmds$points[,2]
ac<- as.data.frame(cbind(a1, a2))

a.nmds.plot<- ggplot()
a.nmds.plot<- a.nmds.plot + geom_vline(x=0,colour="grey50")
a.nmds.plot<- a.nmds.plot + geom_hline(y=0,colour="grey50")
a.nmds.plot<- a.nmds.plot + geom_text(data = ac, aes(x=a1, y=a2, label = a$Approx_year2, colour = a$Lake), size = 7)
a.nmds.plot<- a.nmds.plot + scale_colour_manual(values = c("DAR" = "darkblue", "KB" = "red2"), "Lake")
a.nmds.plot<- a.nmds.plot + theme_bw()
a.nmds.plot<- a.nmds.plot + labs(x= "NMDS 1", y= "NMDS 2")
a.nmds.plot<- a.nmds.plot + theme(axis.text.x = element_text(colour="black", size=16))
a.nmds.plot<- a.nmds.plot + theme(axis.text.y = element_text(colour="black", size=16))
a.nmds.plot<- a.nmds.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
a.nmds.plot<- a.nmds.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
a.nmds.plot<- a.nmds.plot + ggtitle("NMDS - OTUs (Reads (rarefied) non-transformed)")
#Add if want vectors
a.nmds.plot<- a.nmds.plot + geom_path(data = ac, aes(x = a1, y= a2, colour= a$Lake), size =0.4, arrow=arrow(length= unit(0.4,"cm"), ends="first", type="open" ))


##Hierarchical clustering tree - sample by OTU (No. reads)
#Raw data
raw.bray<-vegdist(raw[,5:299], method="bray")

raw.clust<-hclust(raw.bray) 
plot(raw.clust, hang=-1)

raw.clust2<-reorder.hclust(raw.clust, raw.bray)
plot(raw.clust2, hang=-1)

#Using ggdendro
ggdendrogram(raw.clust2, rotate = FALSE, size = 2)

raw.clust3 <- as.dendrogram(raw.clust2)

clust.dat <- dendro_data(raw.clust3, type = "rectangle")

clust.plot <- ggplot() + geom_segment(data = clust.dat$segments, aes(x = x, y = y, xend = xend, yend = yend)) 
clust.plot<- clust.plot + theme_bw()
clust.plot<- clust.plot + geom_text(data = clust.dat$labels, aes(x = x, y = y, label = raw$Approx_year2, colour = raw$Lake), size = 3, vjust = 2)
clust.plot<- clust.plot + scale_colour_manual(values = c("DAR" = "darkblue", "KB" = "red2"), "Lake")


#Rarefied data 
a.bray<-vegdist(a[,5:145], method="bray")

a.clust<-hclust(a.bray) 
plot(a.clust, hang=-1)

a.clust2<-reorder.hclust(a.clust, a.bray)
plot(a.clust2, hang=-1)

#Using ggdendro
ggdendrogram(a.clust2, rotate = FALSE, size = 2)

a.clust3 <- as.dendrogram(a.clust2)

clust.dat <- dendro_data(a.clust3, type = "rectangle")

clust.plot <- ggplot() + geom_segment(data = clust.dat$segments, aes(x = x, y = y, xend = xend, yend = yend)) 
clust.plot<- clust.plot + theme_bw()
clust.plot<- clust.plot + geom_text(data = clust.dat$labels, aes(x = x, y = y, label = a$Approx_year2, colour = a$Lake), size = 3, vjust = 2)
clust.plot<- clust.plot + scale_colour_manual(values = c("DAR" = "darkblue", "KB" = "red2"), "Lake")


##Superimpose clustering objects onto NMDS
#Raw data - 'raw'
raw.bray<-vegdist(raw[,5:299], method="bray")
raw.clust<-hclust(raw.bray)
plot(raw.clust)
raw.grp<- cutree(raw.clust, 2) #
#raw.grp<- cutree(raw.clust, 3) #still only splits out DAR7. 

raw.nmds<-metaMDS(raw[,5:299],distance="bray") 

col<- c("red2", "mediumblue")
col[raw.grp]

plot(raw.nmds, type = "n", display = "sites")
points(raw.nmds, col = col[raw.grp], bg = col[raw.grp], pch = 21)
legend("topright", legend = paste("Cluster", 1:2),
       col = col, pt.bg = col, bty = "n", pch = 21)

plot(raw.nmds, type = "n", display = "sites")
text(raw.nmds, col = col[raw.grp], bg = col[raw.grp], label=raw$Sample_ID)
legend("topright", legend = paste("Cluster", 1:2),
       col = col, pt.bg = col, bty = "n", pch = 21)

#Try with ggplot so can better visualize 
raw.bray<-vegdist(raw[,5:299], method="bray") #distance matrix 
raw.clust<-hclust(raw.bray) #cluster
raw.grp<- as.data.frame(cutree(raw.clust, 2)) #define groups 
raw.nmds<-metaMDS(raw[,5:299],distance="bray") 
raw1 <- raw.nmds$points[,1]
raw2 <- raw.nmds$points[,2]
rawc<- as.data.frame(cbind(raw1, raw2, raw.grp))
colnames(rawc)[3]<- 'grp'
grp<- as.factor(rawc$grp)

raw.nmds.plot<- ggplot()
raw.nmds.plot<- raw.nmds.plot + geom_vline(x=0,colour="grey50")
raw.nmds.plot<- raw.nmds.plot + geom_hline(y=0,colour="grey50")
raw.nmds.plot<- raw.nmds.plot + geom_text(data = rawc, aes(x=raw1, y=raw2, label = raw$Approx_year2, colour = grp), size = 7)
raw.nmds.plot<- raw.nmds.plot + theme_bw()
raw.nmds.plot<- raw.nmds.plot + labs(x= "NMDS 1", y= "NMDS 2")
raw.nmds.plot<- raw.nmds.plot + theme(axis.text.x = element_text(colour="black", size=16))
raw.nmds.plot<- raw.nmds.plot + theme(axis.text.y = element_text(colour="black", size=16))
raw.nmds.plot<- raw.nmds.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
raw.nmds.plot<- raw.nmds.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
raw.nmds.plot<- raw.nmds.plot + ggtitle("NMDS - OTUs (Reads (rarefied) non-transformed, with clusters)")
#Add if want vectors
raw.nmds.plot<- raw.nmds.plot + geom_path(data = rawc, aes(x = raw1, y= raw2, line = raw$Lake), size =0.4, arrow=arrow(length= unit(0.4,"cm"), ends="first", type="open" ))


#Rarefied data - 'a'
a.bray<-vegdist(a[,5:145], method="bray")
a.clust<-hclust(a.bray)
plot(a.clust)
a.grp<- cutree(a.clust, 2) #cutting into 2 groups- if try 3, DAR7 gets lumped on own. 

a.nmds<-metaMDS(a[,5:145],distance="bray") 

col<- c("red2", "mediumblue")
col[a.grp]

plot(a.nmds, type = "n", display = "sites")
points(a.nmds, col = col[a.grp], bg = col[a.grp], pch = 21)
legend("topright", legend = paste("Cluster", 1:2),
       col = col, pt.bg = col, bty = "n", pch = 21)

plot(a.nmds, type = "n", display = "sites")
text(a.nmds, col = col[a.grp], bg = col[a.grp], label=a$Sample_ID)
legend("topright", legend = paste("Cluster", 1:2),
       col = col, pt.bg = col, bty = "n", pch = 21)

#Try with ggplot so can better visualize 
a.bray<-vegdist(a[,5:145], method="bray") #distance matrix 
a.clust<-hclust(a.bray) #cluster
a.grp<- as.data.frame(cutree(a.clust, 2)) #define groups 
a.nmds<-metaMDS(a[,5:145],distance="bray") 
a1 <- a.nmds$points[,1]
a2 <- a.nmds$points[,2]
ac<- as.data.frame(cbind(a1, a2, a.grp))
colnames(ac)[3]<- 'grp'
grp<- as.factor(ac$grp)

a.nmds.plot<- ggplot()
a.nmds.plot<- a.nmds.plot + geom_vline(x=0,colour="grey50")
a.nmds.plot<- a.nmds.plot + geom_hline(y=0,colour="grey50")
a.nmds.plot<- a.nmds.plot + geom_text(data = ac, aes(x=a1, y=a2, label = a$Approx_year2, colour = grp), size = 7)
a.nmds.plot<- a.nmds.plot + theme_bw()
a.nmds.plot<- a.nmds.plot + labs(x= "NMDS 1", y= "NMDS 2")
a.nmds.plot<- a.nmds.plot + theme(axis.text.x = element_text(colour="black", size=16))
a.nmds.plot<- a.nmds.plot + theme(axis.text.y = element_text(colour="black", size=16))
a.nmds.plot<- a.nmds.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
a.nmds.plot<- a.nmds.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
a.nmds.plot<- a.nmds.plot + ggtitle("NMDS - OTUs (Reads (rarefied) non-transformed, with clusters)")
#Add if want vectors
a.nmds.plot<- a.nmds.plot + geom_path(data = ac, aes(x = a1, y= a2, line = a$Lake), size =0.4, arrow=arrow(length= unit(0.4,"cm"), ends="first", type="open" ))


#################MULTIVARIATE REGRESSION TREES##################
##Find breaks in OTU matrix in terms of explanatory variables

#Response matrices - should be read in, transformed and bound to 
#extra lake information the same way as for "Ordination analyses"
raw #sample by OTUs (# of reads, not rarefied) [,5:299]
a #sample by OTUs (# of reads, rarefied) [, 5:145]
#From above in 'ordination analyses' 

#Explanatory variables - concern about too many for n = (check)
raw$Approx_year2  #year (same for a)

#Metal_EFs - extrapolated from geochemical data (from schefferville_lmerdata_MAR2015.xlsx)
#met<- read.table("euk_metalEF.csv",h=T, sep=",")
#raw$met<- met$Metal_EF
#a$met<- met$Metal_EF
#BUT NEED TO FIX THIS- met file will have more rows than 'raw' and 'a'. 

#Use: $Approx_year2, $met (when set-up)


##Multivariate regression trees 

#mvpart - sample by OTUs (# of reads) ('raw' matrix) - BOTH LAKES 
raw.y <- as.matrix(raw[,5:299])
raw.mrt<- mvpart(raw.y ~ raw$Approx_year2, data = raw, cp=0, xv="pick", margin=0.08, xvmult=100, which=4) 
raw.mrt
printcp(raw.mrt) 
plotcp(raw.mrt,minline=TRUE,upper=c("size"))
summary(raw.mrt)
rsq.rpart(raw.mrt)
raw.mrt$cptable
1-raw.mrt$cptable[5,3] #0.42
#Plot
tree.raw.mrt <- MRT(raw.mrt,percent=10,species=colnames(raw.y))
plot(tree.raw.mrt)
summary(tree.raw.mrt)
#Split at 1977 as opposed to the 1955 split visible with eukaryotes. 

#mvpart - sample by OTUs (# of reads) ('a' matrix - rarefied) - BOTH LAKES 
a.y <- as.matrix(a[,5:145])
a.mrt<- mvpart(a.y ~ a$Approx_year2, data = a, cp=0, xv="pick", margin=0.08, xvmult=100, which=4) 
a.mrt
printcp(a.mrt) 
plotcp(a.mrt,minline=TRUE,upper=c("size"))
summary(a.mrt)
rsq.rpart(a.mrt)
a.mrt$cptable
1-a.mrt$cptable[4,3] #0.3
#Plot
tree.a.mrt <- MRT(a.mrt,percent=10,species=colnames(a.y))
plot(tree.a.mrt)
summary(tree.a.mrt)
#Split now at 1920, but this could be because of the samples removed. 

#SPLITTING THE LAKES 
#DAR
raw.DAR<- as.data.frame(subset(raw, Lake == "DAR", drop=T))
a.DAR<- as.data.frame(subset(a, Lake == "DAR", drop=T))

#DAR - raw
raw.DAR.y <- as.matrix(raw.DAR[,5:299])
raw.DAR.mrt<- mvpart(raw.DAR.y ~ raw.DAR$Approx_year2, data = raw.DAR, cp=0, xv="pick", margin=0.08, xvmult=100, which=4) 
raw.DAR.mrt
printcp(raw.DAR.mrt) 
plotcp(raw.DAR.mrt,minline=TRUE,upper=c("size"))
summary(raw.DAR.mrt)
rsq.rpart(raw.DAR.mrt)
raw.DAR.mrt$cptable
1-raw.DAR.mrt$cptable[3,3] #0.52
#Plot
tree.rawDAR.mrt <- MRT(raw.DAR.mrt,percent=10,species=colnames(raw.DAR.y))
plot(tree.rawDAR.mrt)
summary(tree.rawDAR.mrt)
#Split at 1979

#DAR - rarefied
a.DAR.y <- as.matrix(a.DAR[,5:145])
a.DAR.mrt<- mvpart(a.DAR.y ~ a.DAR$Approx_year2, data = a.DAR, cp=0, xv="pick", margin=0.08, xvmult=100, which=4) 
a.DAR.mrt
printcp(a.DAR.mrt) 
plotcp(a.DAR.mrt,minline=TRUE,upper=c("size"))
summary(a.DAR.mrt)
rsq.rpart(a.DAR.mrt)
a.DAR.mrt$cptable
1-a.DAR.mrt$cptable[3,3] #0.39
#Plot
tree.aDAR.mrt <- MRT(a.DAR.mrt,percent=10,species=colnames(a.DAR.y))
plot(tree.aDAR.mrt)
summary(tree.aDAR.mrt)
#Split at 1929. 
  
#KB
raw.KB<- as.data.frame(subset(raw, Lake == "KB", drop=T))
a.KB<- as.data.frame(subset(a, Lake == "KB", drop=T))

#KB - raw
raw.KB.y <- as.matrix(raw.KB[,5:299])
raw.KB.mrt<- mvpart(raw.KB.y ~ raw.KB$Approx_year2, data = raw.KB, cp=0, xv="pick", margin=0.08, xvmult=100, which=4) 
raw.KB.mrt
printcp(raw.KB.mrt) 
plotcp(raw.KB.mrt,minline=TRUE,upper=c("size"))
summary(raw.KB.mrt)
rsq.rpart(raw.KB.mrt)
raw.KB.mrt$cptable
1-raw.KB.mrt$cptable[2,3] #0.46
#Plot
tree.rawKB.mrt <- MRT(raw.KB.mrt,percent=10,species=colnames(raw.KB.y))
plot(tree.rawKB.mrt)
summary(tree.rawKB.mrt)
#No plot, split at 1955. 

#KB - rarefied 
a.KB.y <- as.matrix(a.KB[,5:145])
a.KB.mrt<- mvpart(a.KB.y ~ a.KB$Approx_year2, data = a.KB, cp=0, xv="pick", margin=0.08, xvmult=100, which=4) 
a.KB.mrt
printcp(a.KB.mrt) 
plotcp(a.KB.mrt,minline=TRUE,upper=c("size"))
summary(a.KB.mrt)
rsq.rpart(a.KB.mrt)
a.KB.mrt$cptable
1-a.KB.mrt$cptable[2,3] #0.28
#Plot
tree.aKB.mrt <- MRT(a.KB.mrt,percent=10,species=colnames(a.KB.y))
plot(tree.aKB.mrt)
summary(tree.aKB.mrt)
#No plot, split at 1955. 


#################LINEAR MIXED EFFECT MODELS##################
##Linear-MEMs not useful for eukaryotes, don't
#work on for diatoms for now. 

#################ANOSIM##################
##Statistically test whether there is a signifcant difference
#between groups (vegan pkg) - NB: verified that same results
#are obtained by PAST software. 

#Could do by lake and hclust groups. Time period not important
#for eukaryotes, will need to add hclust groups in. 


#################TEMPORAL BETA-DIVERSITY##################
##Use Legendre functions to compute temporal beta diversity
#between time periods and years. 
#WORK WITH RAW DATA - otherwise run into issues with rarefied data. 

##decompose.D2 function
#####
decompose.D2 <- function(Y1, Y2, den.type=2)
  # Compare two surveys:
  # Decompose the Ruzicka and percentage difference dissimilarities into A, B and C.
  #
  # Parameters --
  # 
  # Y1 : First survey data with sites in rows and species in columns. 
  # Y2 : Second survey data with sites in rows and species in columns. 
  # The sites and species must be the same in the two tables, and in the same order.
  # The files may contain species presence-absence or quantitative abundance data.
  # 
  # den.type -- Denominator type for the indices
#     1 : (A+B+C) as in the Ruzicka dissimilarity
#     2 : (2A+B+C) as in the Percentage difference dissimilarity (alias Bray-Curtis)
#
# Value (output of the function) --
#
# mat1 : A, B and C results, with sites in rows and A, B and C in columns.
# mat2 : A,B,C,D divided by a denominator [either (A+B+C) or (2A+B+C)]; D = (B+C).
#
# Details: Numerical results in output matrices --
# A -- aj is the part of the abundance of species j that is common to the two survey vectors: aj = min(y1j, y2j). A is the sum of the aj values for all species in the functional group under study.
# B -- bj is the part of the abundance of species j that is higher in survey 1 than in survey 2: bj = y1j – y2j. B is the sum of the bj values for all species in the functional group under study.
# C -- cj is the part of the abundance of species j that is higher in survey 2 than in survey 1: cj = y2j – y1j. C is the sum of the cj values for all species in the functional group under study.
#
# Example --
# test1 = matrix(runif(50,0,100),10,5)
# test2 = matrix(runif(50,0,100),10,5)
# (res = decompose.D2(test1, test2, den.type=1))# License: GPL-2
#
# License: GPL-2
# Author:: Pierre Legendre
{
  ### Internal function
  den <- function(A,B,C,den.type) if(den.type==1) den=(A+B+C) else den=(2*A+B+C) #specifying the demoninator to use
  ### End internal function
  #
  n = nrow(Y1)
  p = ncol(Y1)
  if(nrow(Y2)!=n) stop("The data tables do not have the same number of rows") #warning messages in case matrices don't match
  if(ncol(Y2)!=p) stop("The data tables do not have the same number of columns")
  #
  ABC = c("A","B","C")        # A = similarity #column headings for output matrix 1
  ABCD = c("A","B","C","D")   # D = dissimilarity #column headings for output matrix 2
  #
  mat1 = matrix(NA,n,3) #output matrix 1
  colnames(mat1) = ABC
  mat2 = matrix(NA,n,4) #output matrix 2 
  colnames(mat2) = ABCD
  if(!is.null(rownames(Y1))) { 
    rownames(mat1)=rownames(mat2)=rownames(Y1) 
  } else {
    rownames(mat1)=rownames(mat2)=paste("Site",1:n,sep=".") }
  #
  for(i in 1:n) {
    YY = rbind(Y1[i,], Y2[i,])
    A = sum(apply(YY,2,min))
    tmp = YY[1,] - YY[2,]
    B = sum(tmp[tmp>0])
    C = -sum(tmp[tmp<0])
    D = B+C
    mat1[i,] = c(A,B,C)
    mat2[i,] = c(A,B,C,D)/den(A,B,C,den.type) #uses the internal function from above
  }
  #
  list(mat1=mat1, mat2=mat2)
}


#Use 'raw' matrix from ordination and MRT analyses. 
raw.DAR
raw.KB

#DAR
#Make each row a dataframe
#T0 - DAR9 (oldest) - Pre-1850
dar.intT0<- as.data.frame(subset(raw.DAR, Sample_ID == "DAR9", drop=T))

#T1 - DAR8 - 1865
dar.intT1<- as.data.frame(subset(raw.DAR, Sample_ID == "DAR8", drop=T))

#T2 - DAR7 - 1911
dar.intT2<- as.data.frame(subset(raw.DAR, Sample_ID == "DAR7", drop=T))

#T3 - DAR6 - 1947
dar.intT3<- as.data.frame(subset(raw.DAR, Sample_ID == "DAR6", drop=T))

#T4 - DAR5 - 1963
dar.intT4<- as.data.frame(subset(raw.DAR, Sample_ID == "DAR5", drop=T))

#T5 - DAR4 - 1972
dar.intT5<- as.data.frame(subset(raw.DAR, Sample_ID == "DAR4", drop=T))

#T6 - DAR3 - 1986
dar.intT6<- as.data.frame(subset(raw.DAR, Sample_ID == "DAR3", drop=T))

#T7 - DAR2 - 1998
dar.intT7<- as.data.frame(subset(raw.DAR, Sample_ID == "DAR2", drop=T))

#T8 - DAR1 (2012)
dar.intT8<- as.data.frame(subset(raw.DAR, Sample_ID == "DAR1", drop=T))


#Comparisons
#T0-T1 - Pre-1850 to 1865
bdtemp.darT0T1.int<- decompose.D2(dar.intT0[,5:299], dar.intT1[,5:299], den.type=2)  
bdtemp.darT0T1.int.mat2<- as.data.frame(bdtemp.darT0T1.int$mat2) #making result matrix into a dataframe

#T1-T2 - 1865 to 1911
bdtemp.darT1T2.int<- decompose.D2(dar.intT1[,5:299], dar.intT2[,5:299], den.type=2)  
bdtemp.darT1T2.int.mat2<- as.data.frame(bdtemp.darT1T2.int$mat2) 

#T2-T3 - 1911 to 1947
bdtemp.darT2T3.int<- decompose.D2(dar.intT2[,5:299], dar.intT3[,5:299], den.type=2)  
bdtemp.darT2T3.int.mat2<- as.data.frame(bdtemp.darT2T3.int$mat2) 

#T3-T4 - 1947 to 1963
bdtemp.darT3T4.int<- decompose.D2(dar.intT3[,5:299], dar.intT4[,5:299], den.type=2)  
bdtemp.darT3T4.int.mat2<- as.data.frame(bdtemp.darT3T4.int$mat2) 

#T4-T5 - 1963 to 1972
bdtemp.darT4T5.int<- decompose.D2(dar.intT4[,5:299], dar.intT5[,5:299], den.type=2)  
bdtemp.darT4T5.int.mat2<- as.data.frame(bdtemp.darT4T5.int$mat2) 

#T5-T6 - 1972 to 1986
bdtemp.darT5T6.int<- decompose.D2(dar.intT5[,5:299], dar.intT6[,5:299], den.type=2)  
bdtemp.darT5T6.int.mat2<- as.data.frame(bdtemp.darT5T6.int$mat2) 

#T6-T7 - 1986 to 1998
bdtemp.darT6T7.int<- decompose.D2(dar.intT6[,5:299], dar.intT7[,5:299], den.type=2)  
bdtemp.darT6T7.int.mat2<- as.data.frame(bdtemp.darT6T7.int$mat2) 

#T7-T8 - 1998 to 2012
bdtemp.darT7T8.int<- decompose.D2(dar.intT7[,5:299], dar.intT8[,5:299], den.type=2)  
bdtemp.darT7T8.int.mat2<- as.data.frame(bdtemp.darT7T8.int$mat2) 

#KB
#Make each row a dataframe
#T0 - KB6 (oldest) - Pre-1850
kb.intT0<- as.data.frame(subset(raw.KB, Sample_ID == "KB6", drop=T))

#T1 - KB5 - Pre-1850
kb.intT1<- as.data.frame(subset(raw.KB, Sample_ID == "KB5", drop=T))

#T2 - KB4 - 1866
kb.intT2<- as.data.frame(subset(raw.KB, Sample_ID == "KB4", drop=T))

#T3 - KB3 - 1928
kb.intT3<- as.data.frame(subset(raw.KB, Sample_ID == "KB3", drop=T))

#T4 - KB2 - 1982
kb.intT4<- as.data.frame(subset(raw.KB, Sample_ID == "KB2", drop=T))

#T5 - KB1 (2012)
kb.intT5<- as.data.frame(subset(raw.KB, Sample_ID == "KB1", drop=T))

#Comparisons
#T0-T1- Pre1850(1) to Pre1850(2)
bdtemp.kbT0T1.int<- decompose.D2(kb.intT0[,5:299], kb.intT1[,5:299], den.type=2)  
bdtemp.kbT0T1.int.mat2<- as.data.frame(bdtemp.kbT0T1.int$mat2) 

#T1-T2- Pre1850(2) to 1866
bdtemp.kbT1T2.int<- decompose.D2(kb.intT1[,5:299], kb.intT2[,5:299], den.type=2)  
bdtemp.kbT1T2.int.mat2<- as.data.frame(bdtemp.kbT1T2.int$mat2) 

#T2-T3- 1866 ro 1928
bdtemp.kbT2T3.int<- decompose.D2(kb.intT2[,5:299], kb.intT3[,5:299], den.type=2)  
bdtemp.kbT2T3.int.mat2<- as.data.frame(bdtemp.kbT2T3.int$mat2) 

#T3-T4- 1928 to 1982
bdtemp.kbT3T4.int<- decompose.D2(kb.intT3[,5:299], kb.intT4[,5:299], den.type=2)  
bdtemp.kbT3T4.int.mat2<- as.data.frame(bdtemp.kbT3T4.int$mat2) 

#T4-T5- 1982 to 2012
bdtemp.kbT4T5.int<- decompose.D2(kb.intT4[,5:299], kb.intT5[,5:299], den.type=2)  
bdtemp.kbT4T5.int.mat2<- as.data.frame(bdtemp.kbT4T5.int$mat2) 


#Explore results with plots
beta<- read.csv("diat_temporalBD_15sRAWr.csv") #Remember, with RAW data. 
beta.long<- melt(beta, id.vars=c("Lake", "Comparison_type", "Comparison", "Comparison_years", "Comparison_order"))
colnames(beta.long)[6]<- 'Beta_component'
colnames(beta.long)[7]<- 'Prop'
as.factor(beta.long$Comparison_order)


#DAR
beta.long.DAR<- as.data.frame(subset(beta.long, Lake == "DAR", drop=T))
beta.long.DAR<- as.data.frame(subset(beta.long.DAR, Beta_component == "Total_beta" | Beta_component == "Species_loss" | Beta_component == "Species_gain", drop=T))
#Don't need Prop_loss and Prop_gain columns 

betaDAR.plot<- ggplot(beta.long.DAR, aes(x=Comparison_order, y=Prop, colour = Beta_component)) + geom_point(size=4) + geom_path(size=2)
betaDAR.plot<- betaDAR.plot + theme_bw() + labs(x= "Comparison", y= "Beta")
betaDAR.plot<- betaDAR.plot + theme(axis.text.x = element_text(colour="black", size=16, angle=45))
betaDAR.plot<- betaDAR.plot + theme(axis.text.y = element_text(colour="black", size=16))
betaDAR.plot<- betaDAR.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
betaDAR.plot<- betaDAR.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
betaDAR.plot<- betaDAR.plot + scale_x_discrete(labels=c("Pre1850-1865","1865-1911","1911-1947","1947-1963", "1963-1972", "1972-1986", "1986-1998", "1998-2012"))
betaDAR.plot<- betaDAR.plot + ggtitle("Dauriat")


#KB
beta.long.KB<- as.data.frame(subset(beta.long, Lake == "KB", drop=T))
beta.long.KB<- as.data.frame(subset(beta.long.KB, Beta_component == "Total_beta" | Beta_component == "Species_loss" | Beta_component == "Species_gain", drop=T))

betaKB.plot<- ggplot(beta.long.KB, aes(x=Comparison_order, y=Prop, colour = Beta_component)) + geom_point(size=4) + geom_path(size=2)
betaKB.plot<- betaKB.plot + theme_bw() + labs(x= "Comparison", y= "Beta")
betaKB.plot<- betaKB.plot + theme(axis.text.x = element_text(colour="black", size=16, angle=45))
betaKB.plot<- betaKB.plot + theme(axis.text.y = element_text(colour="black", size=16))
betaKB.plot<- betaKB.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
betaKB.plot<- betaKB.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
betaKB.plot<- betaKB.plot + scale_x_discrete(labels=c("Pre1850-Pre1850","Pre1850-1866","1866-1928","1928-1982", "1982-2012"))
betaKB.plot<- betaKB.plot + ggtitle("Knob")

