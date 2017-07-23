##################################################################################################
####SCHEFFERVILLE NGS: Sept 2015 work                                                            #
##################################################################################################

##In this script: Preliminary analyses of cleaned eukaryote data from
#Sept 2015 PANAM output. 
#Previous script: Script_ecological_inferences.R 
#R version: 3.1.2 (Pumpkin Helmet)

##Last update: Oct. 28, 2015
##Associated workspace: 
##Associated markdown: 
##Associated .txt of R script: 
##Github: 

##################################################################################################

##Working with eukaryote data that has already been matched up with 
#taxonomy and family-/functional group level matches in ecological 
#inferences script. Singletons removed at outset. 

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

library(devtools) #compiler for archived packages #need to DL each time?
#Also installed Rtools (but not a package that needs to be loaded)

#library(rdaTEST) #not available for R. 3.1.2 

#mvpart: install.packages("C:/Users/Winegardner/Downloads/mvpart_1.6-2.tar.gz", repos = NULL, type = "source")
#MVPARTwrap: install.packages("C:/Users/Winegardner/Downloads/MVPARTwrap_0.1-9.tar.gz", repos = NULL, type = "source")

#################NORMY FUNCTION##################
##Function to "normalize" data by the minimum number of reads (rarefied read #).  
#Remove sample 10 (the last DAR sample as affecting the normalization, too low)

a <- read.table("euk_nosing_ss_EMBL.csv",h=T, sep=",") #eukaryote data, no singletons, matched to taxonomy
#Remove column 11 (Sed10)
a<- a[,-11]
data.frame(colnames(a))
min(colSums(a[,2:16])) #for minimum number of reads, 741 if remove sample 10

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

b <- normy(data=a,number_of_samples=15,number_of_reads=741)
write.table(b,file="euk_nosing_15s741r.txt",sep="\t",row.names=F) #labelled with 15 samples and 741 reads.

#################RAREFACTION CURVES##################
##Create rarefaction curves to see aysmptotes for rarefied read counts. 

a<-read.table("euk_nosing_15s741r.txt",h=T,row.names=1) #uses output from normy
a1<-as.data.frame(t(a[,1:15]))
S<-specnumber(a1)
raremax<-min(rowSums(a1))
Srare<-rarefy(a1,raremax)
plot(S,Srare,xlab="ObservedNo.ofSpecies",ylab="RarefiedNo.ofSpecies")
abline(0,1)
rarecurve(a1,step=100,sample=592,col="blue",cew=0.6)

#################METRY FUNCTION:CALCULATING RICHNESS AND DIVERSITY METRICS##################
##Function to calculate variety of diversity metrics on read data. 

a <- read.table("euk_nosing_15s741r.txt",h=T) #use rarefied data (output from normy) 
#could also produce using non-rarefied data. 
data.frame(colnames(a))

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

b<-metry(data=a,number_of_samples=15)
write.table(b[,1:11],file="euk_nosing_15s741r_metrics.txt",sep="\t",row.names=F)
#only need 1st 11 columns. 

#Re-fill in "expanded" version with dates if change. 

#################TAXY:TAXO GROUPING##################
##Group reads and OTUs into functional groups (as defined by either EMBL2 or EMBL3).
#Produce matrices of # of reads and # of OTUs for each of these functional groupings.

a<-read.table("euk_nosing_15s741r.txt",h=T) #rarefied data (output from normy).
data.frame(colnames(a))

#Produce table with # of reads for each functional group, for each interval. 
#using taxonomic level = EMBL3
taxy_nbreads = function(data, number_of_samples, column_of_taxo_level, name_of_taxo_level){
  
  data1 <- subset(data, data[,length(data)] !="0")    
  data2 <- aggregate(data1[,2:(number_of_samples+1)], by = list(typ2 = data1[,column_of_taxo_level]), sum) 
  data3 <- aggregate(data1[,2], by = list(typ2 = data1[,column_of_taxo_level]), length) 
  abund_taxo <- colSums(t(data2[,2:(number_of_samples+1)])) 
  data4 <- cbind(data2,data3, abund_taxo)
  data5 <- data4[,-(number_of_samples+2)]
  colnames(data5) <- c(name_of_taxo_level, colnames(data5[,2:(number_of_samples+1)]), "Nb_OTUs", "Abund")
  data5}
b <- taxy_nbreads(data=a,number_of_samples=15,column_of_taxo_level=19,name_of_taxo_level="EMBL3")
write.table(b,file="euk_nosing_15s741r_EMBL3_nbreads.txt",sep="\t",row.names=F)

#Produce table with # of reads for each functional group, for each interval. 
#using taxonomic level = EMBL 2
b <- taxy_nbreads(data=a,number_of_samples=15,column_of_taxo_level=20,name_of_taxo_level="EMBL2")
write.table(b,file="euk_nosing_15s741r_EMBL2_nbreads.txt",sep="\t",row.names=F)


#Produce table with # of OTUs for each functional group, for each interval. 
#using taxonomic level = EMBL3
taxy_nbotus = function(data, number_of_samples, column_of_taxo_level, name_of_taxo_level){
  data1 <- subset(data, data[,length(data)] !="0")    
  data2 <- list()  
  for(i in 2:(number_of_samples+1))
  {data2[[i]] <- as.numeric(data1[,i] > 0)}
  data3 <- matrix(unlist(data2), ncol = number_of_samples)
  data_bin <- cbind(data1[,1], data3, data1[,(number_of_samples+2):length(data1)])
  colnames(data_bin) <- colnames(data)
  data_bin2 <- aggregate(data_bin[,2:(number_of_samples+1)], by = list(typ2 = data_bin[,column_of_taxo_level]), sum) 
  data_bin3 <- aggregate(data_bin[,2], by = list(typ2 = data_bin[,column_of_taxo_level]), length) 
  data4 <- cbind(data_bin2,data_bin3)
  data5 <- data4[,-(number_of_samples+2)]
  colnames(data5) <- c(name_of_taxo_level, colnames(data1[,2:(number_of_samples+1)]), "Total_Nb_OTUs")
  data5}
c <- taxy_nbotus(data=a,number_of_samples=15,column_of_taxo_level=19,name_of_taxo_level="EMBL3")
write.table(c,file="euk_nosing_15s741r_EMBL3_nbotus.txt",sep="\t",row.names=F)

#Produce table with # of OTUs for each functional group, for each interval. 
#using taxonomic level = EMBL2
c <- taxy_nbotus(data=a,number_of_samples=15,column_of_taxo_level=20,name_of_taxo_level="EMBL2")
write.table(c,file="euk_nosing_15s741r_EMBL2_nbotus.txt",sep="\t",row.names=F)

#################HISTOGRAMS EMBL3 and EMBL2##################
##Produce some exploratory histograms using functional groupings EMBL3 and 2.

##EMBL3
#Run histogram first with one, then with other. 
a<-read.table("euk_nosing_15s741r_EMBL3_nbreads.txt",h=T) #when run with # of reads, normalized to 741
a<-read.table("euk_nosing_15s741r_EMBL3_nbotus.txt",h=T) #when run, total OTUs will be cap (not all out of same total)

a1<-sub("Apicomplexa","darkgoldenrod1",a[,1])
a1<-sub("Ciliophora","darkblue",a1)
a1<-sub("unclassified_Ciliophora","darkblue",a1)
a1<-sub("Colpodellidae","darkcyan",a1)
a1<-sub("Dinophyceae","green3",a1)
a1<-sub("ES_Alveolata","royalblue",a1)
a1<-sub("Perkinsea","lightblue4",a1)
a1<-sub("Voromonas","blue4",a1)
a1<-sub("unclassified_Alveolata","royalblue",a1)
a1<-sub("Discosea","indianred4",a1)
a1<-sub("Mycetozoa","indianred4",a1)
a1<-sub("Tubulinea","indianred4",a1)
a1<-sub("unclassified_Amoebozoa","indianred4",a1)
a1<-sub("Cryptomonadales","yellow",a1)
a1<-sub("Cryptophyta_2","yellow",a1)
a1<-sub("Cryptophyta_3","yellow",a1)
a1<-sub("Cryptophyta_4","yellow",a1)
a1<-sub("Cryptophyta_5","yellow",a1)
a1<-sub("ES_Cryptophyta","yellow",a1)
a1<-sub("Pyrenomonadales","yellow",a1)
a1<-sub("unclassified_Cryptophyta","yellow",a1)
a1<-sub("Euglenida","mediumpurple2",a1)
a1<-sub("Kinetoplastida","mediumpurple2",a1)
a1<-sub("Diplomonadida","skyblue3",a1)
a1<-sub("Coccolithales","turquoise2",a1)
a1<-sub("Coccosphaerales","turquoise2",a1)
a1<-sub("ES_Haptophyceae","turquoise2",a1)
a1<-sub("Isochrysidales","turquoise2",a1)
a1<-sub("Pavlovales","turquoise2",a1)
a1<-sub("Phaeocystales","turquoise2",a1)
a1<-sub("Prymnesiales","turquoise2",a1)
a1<-sub("Syracosphaerales","turquoise2",a1)
a1<-sub("unclassified_Haptophyceae","turquoise2",a1)
a1<-sub("Zygodiscales","turquoise2",a1)
a1<-sub("Schizopyrenida","gray",a1)
a1<-sub("Nucleariidae","rosybrown1",a1)
a1<-sub("Choanoflagellida","gray",a1)
a1<-sub("Fungi","chocolate4",a1)
a1<-sub("Opisthokonta_incertae_sedis","gray",a1)
a1<-sub("Trichomonadida","darkolivegreen3",a1)
a1<-sub("Honigbergiellida","darkolivegreen3",a1)
a1<-sub("Hypotrichomonadida","darkolivegreen3",a1)
a1<-sub("Cercozoa","seagreen3",a1)
a1<-sub("ES_Rhizaria","purple2",a1)
a1<-sub("Foraminifera","purple2",a1)
a1<-sub("Bangiophyceae","red",a1)
a1<-sub("Florideophyceae","red",a1)
a1<-sub("Bacillariophyta","yellowgreen",a1)
a1<-sub("Bicosoecida","khaki1",a1)
a1<-sub("Bolidophyceae","khaki2",a1)
a1<-sub("Chrysophyceae","seagreen1",a1)
a1<-sub("Dictyochophyceae","orangered",a1)
a1<-sub("ES_Stramenopiles","orangered",a1)
a1<-sub("Eustigmatophyceae","orangered",a1)
a1<-sub("Hyphochytriomycetes","orangered",a1)
a1<-sub("Labyrinthulomycetes","orangered",a1)
a1<-sub("Labyrinthulida","orangered",a1)
a1<-sub("Leukarachnion","orangered",a1)
a1<-sub("Oikomonadaceae","orangered",a1)
a1<-sub("Oomycetes","lightsalmon4",a1)
a1<-sub("Pelagophyceae","orangered",a1)
a1<-sub("Phaeothamniophyceae","orangered",a1)
a1<-sub("Pinguiophyceae","orangered",a1)
a1<-sub("Pirsonia","orangered",a1)
a1<-sub("PX_clade","orangered",a1)
a1<-sub("Raphidophyceae","orangered",a1)
a1<-sub("Slopalinida","orangered",a1)
a1<-sub("Raphidophyceae","orangered",a1)
a1<-sub("MAST-1","orangered",a1)
a1<-sub("MAST-2","orangered",a1)
a1<-sub("MAST-3","orangered",a1)
a1<-sub("MAST-4","orangered",a1)
a1<-sub("MAST-7","orangered",a1)
a1<-sub("MAST-6","orangered",a1)
a1<-sub("MAST-7","orangered",a1)
a1<-sub("MAST-12","orangered",a1)
a1<-sub("Synurophyceae","orangered",a1)
a1<-sub("unclassified_Stramenopiles","orangered",a1)
a1<-sub("Chlorophyta","green4",a1)
a1<-sub("ES_Viridiplantae","green2",a1)
a1<-sub("Streptophyta","lightpink4",a1)
a1<-sub("unclassified_Viridiplantae","green2",a1)
a1<-sub("unclassified_Rhizaria","purple2",a1)
a1

barplot(as.matrix(a[,2:16]),las=2,col=a1, main="No. OTUs - EMBL3") #change title as needed
#Title options:
#"No. Reads - EMBL3"
#"No. OTUs - EMBL3"
#See Eric's legend for colours. 

##Make stratigraphy of EMBL3 in ggplot2

#EMBL3 with # of reads 
a<-read.table("euk_nosing_15s741r_EMBL3_nbreads.txt",h=T)

#Transpose so sediments in long format
a<- a[,1:16] #only need OTUs and samples
a<- melt(a) #make long-format
a<- dcast(a, variable ~ EMBL3) #cast back into transposed matrix. 
colnames(a) [1]<- "Sample_ID"
rownames(a)<- as.character(a[,1])

#Add in columns for lake and time period
div <- read.csv("euk_nosing_15s741r_metrics_expanded.csv") #read in data that has this information
b<- as.data.frame(cbind(div[,2:5], a[,2:37]))
#Melt back into longformat so can use in ggplot. 
c<- melt(b, id.vars = c("Lake", "Sample_ID", "Approx_year", "Approx_year2"))
colnames(c)[5]<- 'EMBL3'
colnames(c)[6]<- 'No.reads'

#Make histogram
d<- ggplot(c, aes(x=Approx_year2, y=No.reads, fill=EMBL3)) + geom_bar(stat = 'identity', width = 5.00) + coord_flip() + facet_wrap(~Lake)
d<- d + scale_fill_manual(values = c("Apicomplexa" = "darkgoldenrod1", "Bacillariophyta" = "darkblue",
                          "Bicosoecida" = "darkslateblue", "Bolidophyceae" = "khaki", "Cercozoa" = "seagreen3", 
                          "Chlorophyta" = "green", "Choanoflagellida" = "orangered3", "Chrysophyceae" = "red1", 
                          "Ciliophora" = "mediumvioletred", "Cryptomonadales" = "yellow2", "Cryptophyta_3" = "yellow2",
                          "Cryptophyta_5" = "yellow2", "Dinophyceae" = "turquoise4", "Discosea" = "wheat4",
                          "ES_Alveolata" = "slategrey", "ES_Haptophyceae" = "slategrey", "ES_Rhizaria" = "slategrey",
                          "ES_Stramenopiles" = "slategrey", "Euglenida" = "olivedrab3", "Eustigmatophyceae" = "limegreen",
                          "Fungi" = "forestgreen", "Hypotrichomonadida" = "deeppink", "Kinetoplastida" = "darksalmon",
                          "Labyrinthulida" = "dimgray", "MAST-2" = "firebrick1", "MAST-3" = "firebrick1", 
                          "Oikomonadaceae" = "chocolate2", "Oomycetes" = "cadetblue4", "Opisthokonta_incertae_sedis" = "blue",
                          "Perkinsea" = "brown", "Prymnesiales" = "azure4", "Pyrenomonadales" = "deeppink4",
                          "Raphidophyceae" = "gray11", "Streptophyta" = "darkturquoise", "Synurophyceae" = "darkslategray1",
                          "unclassified_Stramenopiles" = "darkred"))
d<- d + labs(x= "Year", y= "Number of reads") + theme_bw()
d<- d + theme(axis.text.x = element_text(colour="black", size=16))
d<- d + theme(axis.text.y = element_text(colour="black", size=16))
d<- d + theme(axis.title.x = element_text(size = rel(2), angle=00))
d<- d + theme(axis.title.y = element_text(size = rel(2), angle=90))
#Confirm that KB is normalized to rarefied # of reads. 

#EMBL with # of OTUs
a<-read.table("euk_nosing_15s741r_EMBL3_nbotus.txt",h=T)

#Transpose so sediments in long format
a<- a[,1:16] #only need OTUs and samples
a<- melt(a) #make long-format
a<- dcast(a, variable ~ EMBL3) #cast back into transposed matrix. 
colnames(a) [1]<- "Sample_ID"
rownames(a)<- as.character(a[,1])

#Add in columns for lake and time period
div <- read.csv("euk_nosing_15s741r_metrics_expanded.csv") #read in data that has this information
b<- as.data.frame(cbind(div[,2:5], a[,2:37]))
#Melt back into longformat so can use in ggplot. 
c<- melt(b, id.vars = c("Lake", "Sample_ID", "Approx_year", "Approx_year2"))
colnames(c)[5]<- 'EMBL3'
colnames(c)[6]<- 'No.OTUs'

#Make histogram
d<- ggplot(c, aes(x=Approx_year2, y=No.OTUs, fill=EMBL3)) + geom_bar(stat = 'identity', width = 5.00) + coord_flip() + facet_wrap(~Lake)
d<- d + scale_fill_manual(values = c("Apicomplexa" = "darkgoldenrod1", "Bacillariophyta" = "darkblue",
                                     "Bicosoecida" = "darkslateblue", "Bolidophyceae" = "khaki", "Cercozoa" = "seagreen3", 
                                     "Chlorophyta" = "green", "Choanoflagellida" = "orangered3", "Chrysophyceae" = "red1", 
                                     "Ciliophora" = "mediumvioletred", "Cryptomonadales" = "yellow2", "Cryptophyta_3" = "yellow2",
                                     "Cryptophyta_5" = "yellow2", "Dinophyceae" = "turquoise4", "Discosea" = "wheat4",
                                     "ES_Alveolata" = "slategrey", "ES_Haptophyceae" = "slategrey", "ES_Rhizaria" = "slategrey",
                                     "ES_Stramenopiles" = "slategrey", "Euglenida" = "olivedrab3", "Eustigmatophyceae" = "limegreen",
                                     "Fungi" = "forestgreen", "Hypotrichomonadida" = "deeppink", "Kinetoplastida" = "darksalmon",
                                     "Labyrinthulida" = "dimgray", "MAST-2" = "firebrick1", "MAST-3" = "firebrick1", 
                                     "Oikomonadaceae" = "chocolate2", "Oomycetes" = "cadetblue4", "Opisthokonta_incertae_sedis" = "blue",
                                     "Perkinsea" = "brown", "Prymnesiales" = "azure4", "Pyrenomonadales" = "deeppink4",
                                     "Raphidophyceae" = "gray11", "Streptophyta" = "darkturquoise", "Synurophyceae" = "darkslategray1",
                                     "unclassified_Stramenopiles" = "darkred"))
d<- d + labs(x= "Year", y= "Number of OTUs") + theme_bw()
d<- d + theme(axis.text.x = element_text(colour="black", size=16))
d<- d + theme(axis.text.y = element_text(colour="black", size=16))
d<- d + theme(axis.title.x = element_text(size = rel(2), angle=00))
d<- d + theme(axis.title.y = element_text(size = rel(2), angle=90))


##EMBL2
#Use one or other to run 
a<-read.table("euk_nosing_15s741r_EMBL2_nbreads.txt",h=T) #when run with # of reads, normalized to 741
a<-read.table("euk_nosing_15s741r_EMBL2_nbotus.txt",h=T) #when run, total OTUs will be cap (not all out of same total)

a1<-sub("Alveolata","darkgoldenrod1",a[,1])
a1<-sub("Amoebozoa","darkblue",a1)
a1<-sub("Cryptophyta","olivedrab1",a1)
a1<-sub("Euglenozoa","green3",a1)
a1<-sub("Haptophyceae","royalblue",a1)
a1<-sub("Hexamitidae","lightblue4",a1)
a1<-sub("Opisthokonta","blue4",a1)
a1<-sub("Parabasalia","pink",a1)
a1<-sub("Rhizaria","indianred4",a1)
a1<-sub("Stramenopiles","darkviolet",a1)
a1<-sub("Viridiplantae","green1",a1)
a1

barplot(as.matrix(a[,2:16]),las=2,col=a1, main="No. OTUs - EMBL2") #change title as needed. 
#Title options:
#"No. reads - EMBL2"
#"No. OTUs - EMBL3" 
#See colours above for legend. 

##Make stratigraphy of EMBL2 in ggplot2

#EMBL2 with # of reads 
a<-read.table("euk_nosing_15s741r_EMBL2_nbreads.txt",h=T)

#Transpose so sediments in long format
a<- a[,1:16] #only need OTUs and samples
a<- melt(a) #make long-format
a<- dcast(a, variable ~ EMBL2) #cast back into transposed matrix. 
colnames(a) [1]<- "Sample_ID"
rownames(a)<- as.character(a[,1])

#Add in columns for lake and time period
div <- read.csv("euk_nosing_15s741r_metrics_expanded.csv") #read in data that has this information
b<- as.data.frame(cbind(div[,2:5], a[,2:11]))
#Melt back into longformat so can use in ggplot. 
c<- melt(b, id.vars = c("Lake", "Sample_ID", "Approx_year", "Approx_year2"))
colnames(c)[5]<- 'EMBL2'
colnames(c)[6]<- 'No.reads'

#Make histogram
d<- ggplot(c, aes(x=Approx_year2, y=No.reads, fill=EMBL2)) + geom_bar(stat = 'identity', width = 5.00) + coord_flip() + facet_wrap(~Lake)
d<- d + scale_fill_manual(values = c("Alveolata" = "firebrick3", "Amoebozoa" = "darkslateblue",
                                     "Cryptophyta" = "goldenrod3", "Euglenozoa" = "forestgreen",
                                     "Haptophyceae" = "darkred", "Opisthokonta" = "darkmagenta",
                                     "Parabasalia" = "dodgerblue", "Rhizaria" = "darkorange3",
                                     "Stramenopiles" = "darkolivegreen3", "Viridiplantae" = "lawngreen"))
d<- d + labs(x= "Year", y= "Number of reads") + theme_bw()
d<- d + theme(axis.text.x = element_text(colour="black", size=16))
d<- d + theme(axis.text.y = element_text(colour="black", size=16))
d<- d + theme(axis.title.x = element_text(size = rel(2), angle=00))
d<- d + theme(axis.title.y = element_text(size = rel(2), angle=90))
#Confirm that KB is normalized to rarefied # of reads. 

#EMBL2 with # of OTUs
a<-read.table("euk_nosing_15s741r_EMBL2_nbotus.txt",h=T)

#Transpose so sediments in long format
a<- a[,1:16] #only need OTUs and samples
a<- melt(a) #make long-format
a<- dcast(a, variable ~ EMBL2) #cast back into transposed matrix. 
colnames(a) [1]<- "Sample_ID"
rownames(a)<- as.character(a[,1])

#Add in columns for lake and time period
div <- read.csv("euk_nosing_15s741r_metrics_expanded.csv") #read in data that has this information
b<- as.data.frame(cbind(div[,2:5], a[,2:11]))
#Melt back into longformat so can use in ggplot. 
c<- melt(b, id.vars = c("Lake", "Sample_ID", "Approx_year", "Approx_year2"))
colnames(c)[5]<- 'EMBL2'
colnames(c)[6]<- 'No.OTUs'

#Make histogram
d<- ggplot(c, aes(x=Approx_year2, y=No.OTUs, fill=EMBL2)) + geom_bar(stat = 'identity', width = 5.00) + coord_flip() + facet_wrap(~Lake)
d<- d + scale_fill_manual(values = c("Alveolata" = "firebrick3", "Amoebozoa" = "darkslateblue",
                                     "Cryptophyta" = "goldenrod3", "Euglenozoa" = "forestgreen",
                                     "Haptophyceae" = "darkred", "Opisthokonta" = "darkmagenta",
                                     "Parabasalia" = "dodgerblue", "Rhizaria" = "darkorange3",
                                     "Stramenopiles" = "darkolivegreen3", "Viridiplantae" = "lawngreen"))
d<- d + labs(x= "Year", y= "Number of OTUs") + theme_bw()
d<- d + theme(axis.text.x = element_text(colour="black", size=16))
d<- d + theme(axis.text.y = element_text(colour="black", size=16))
d<- d + theme(axis.title.x = element_text(size = rel(2), angle=00))
d<- d + theme(axis.title.y = element_text(size = rel(2), angle=90))


#################VISUALIZE DIVERSITY METRICS##################
##Load the diversity metrics from metry function and visualize changes
#through cores and time. 

div <- read.csv("euk_nosing_15s741r_metrics_expanded.csv") #Exported metrics with lake info and approx dates added. 

#Subset into lac Dauriat and Knob
div.DAR<- as.data.frame(subset(div, Lake == "DAR", drop=T))
div.KB<- as.data.frame(subset(div, Lake == "KB", drop=T))

#DAR plots
#Substitute y = and ylab =
#Shannon
#Simpson
#Invsimpson
#Pielou
#Chao1
#Chao2
#ACE
#ICE
#Goods_coverage

#Dating years error
error<- aes(xmax = div.DAR$Approx_yearCIC + div.DAR$Est_year_error, xmin = div.DAR$Approx_yearCIC - div.DAR$Est_year_error)

#DAR- Shannon
dar.plot1<- ggplot(div.DAR, aes(x=Approx_yearCIC, y=Shannon)) + geom_path(size=0.75, colour="darkblue")
dar.plot1<- dar.plot1 + geom_point(aes(x=Approx_yearCIC, y=Shannon, shape=factor(Year_measurement)), size=4, colour="darkblue")
dar.plot1<- dar.plot1 + scale_shape_manual(values = c(16,1))
dar.plot1<- dar.plot1 + labs(x= "Estimated year", y= "Shannon diversity") + theme_bw() + theme(legend.position="none")
dar.plot1<- dar.plot1 + geom_errorbarh(error, width=0.25)
dar.plot1<- dar.plot1 + theme(axis.text.x = element_text(colour="black", size=16))
dar.plot1<- dar.plot1 + theme(axis.text.y = element_text(colour="black", size=16))
dar.plot1<- dar.plot1 + theme(axis.title.x = element_text(size = rel(2), angle=00))
dar.plot1<- dar.plot1 + theme(axis.title.y = element_text(size = rel(2), angle=90))
dar.plot1<- dar.plot1 + annotate("rect", xmin=1939, xmax=1977, ymin=1, ymax=5, alpha=0.2) #rectangle highlights mining time period
#rectangle not perfect match for all metrics. 
dar.plot1<- dar.plot1 + annotate("text", x=1920, y=5, label = "(a)", size=12)

#DAR- No. OTUs (rarefied)
dar.plot2<- ggplot(div.DAR, aes(x=Approx_yearCIC, y=Nb_otus_rarefied)) + geom_path(size=0.75, colour="darkblue")
dar.plot2<- dar.plot2 + geom_point(aes(x=Approx_yearCIC, y=Nb_otus_rarefied, shape=factor(Year_measurement)), size=4, colour="darkblue")
dar.plot2<- dar.plot2 + scale_shape_manual(values = c(16,1))
dar.plot2<- dar.plot2 + labs(x= "Estimated year", y= "# OTUs") + theme_bw() + theme(legend.position="none")
dar.plot2<- dar.plot2 + geom_errorbarh(error, width=0.25)
dar.plot2<- dar.plot2 + theme(axis.text.x = element_text(colour="black", size=16))
dar.plot2<- dar.plot2 + theme(axis.text.y = element_text(colour="black", size=16))
dar.plot2<- dar.plot2 + theme(axis.title.x = element_text(size = rel(2), angle=00))
dar.plot2<- dar.plot2 + theme(axis.title.y = element_text(size = rel(2), angle=90))
dar.plot2<- dar.plot2 + annotate("rect", xmin=1939, xmax=1977, ymin=1, ymax=250, alpha=0.2)

#dar.plotcombo<- grid.arrange(dar.plot1, dar.plot2, nrow=1)

#KB plots
#Substitute y = and ylab =
#Shannon
#Simpson
#Invsimpson
#Pielou
#Chao1
#Chao2
#ACE
#ICE
#Goods_coverage

#Dating years error
errorKB<- aes(xmax = div.KB$Approx_yearCIC + div.KB$Est_year_error, xmin = div.KB$Approx_yearCIC - div.KB$Est_year_error)

#KB- Shannon
kb.plot1<- ggplot(div.KB, aes(x=Approx_yearCIC, y=Shannon)) + geom_path(size=0.75, colour="darkred")
kb.plot1<- kb.plot1 + geom_point(aes(x=Approx_yearCIC, y=Shannon, shape=factor(Year_measurement)), size=4, colour="darkred")
kb.plot1<- kb.plot1 + scale_shape_manual(values = c(16,1))
kb.plot1<- kb.plot1 + labs(x= "Estimated year", y= "Shannon diversity") + theme_bw() + theme(legend.position="none")
kb.plot1<- kb.plot1 + geom_errorbarh(errorKB, width=0.25)
kb.plot1<- kb.plot1 + theme(axis.text.x = element_text(colour="black", size=16))
kb.plot1<- kb.plot1 + theme(axis.text.y = element_text(colour="black", size=16))
kb.plot1<- kb.plot1 + theme(axis.title.x = element_text(size = rel(2), angle=00))
kb.plot1<- kb.plot1 + theme(axis.title.y = element_text(size = rel(2), angle=90))
kb.plot1<- kb.plot1 + annotate("rect", xmin=1939, xmax=1977, ymin=1, ymax=5, alpha=0.2) 
#rectangle not perfect match for all metrics.
kb.plot1<- kb.plot1 + annotate("text", x=1736, y=5, label = "(b)", size=12)

#KB- No. OTUs
kb.plot2<- ggplot(div.KB, aes(x=Approx_yearCIC, y=Nb_otus_rarefied)) + geom_path(size=0.75, colour="darkred")
kb.plot2<- kb.plot2 + geom_point(aes(x=Approx_yearCIC, y=Nb_otus_rarefied, shape=factor(Year_measurement)), size=4, colour="darkred")
kb.plot2<- kb.plot2 + scale_shape_manual(values = c(16,1))
kb.plot2<- kb.plot2 + labs(x= "Estimated year", y= "# OTUs") + theme_bw() + theme(legend.position="none")
kb.plot2<- kb.plot2 + geom_errorbarh(errorKB, width=0.25)
kb.plot2<- kb.plot2 + theme(axis.text.x = element_text(colour="black", size=16))
kb.plot2<- kb.plot2 + theme(axis.text.y = element_text(colour="black", size=16))
kb.plot2<- kb.plot2 + theme(axis.title.x = element_text(size = rel(2), angle=00))
kb.plot2<- kb.plot2 + theme(axis.title.y = element_text(size = rel(2), angle=90))
kb.plot2<- kb.plot2 + annotate("rect", xmin=1939, xmax=1977, ymin=1, ymax=263, alpha=0.2) 


#Combo
alpha.plot<- grid.arrange(dar.plot1, dar.plot2, kb.plot1, kb.plot2, nrow=2)


#################ORDINATION ANALYSES##################
##Generate some basic ordinations in order to explore variation in
#community assemblages and structure. 

#Matrices
a<- read.table("euk_nosing_15s741r.txt",h=T) #OTU by sample matrix with rarefied reads
#columns 2:16 are the sediment samples. 

b<- read.table("euk_nosing_15s741r_EMBL3_nbreads.txt",h=T) #EMBL3 functional group by sample matrix with rarefied reads
#columns 2:16 are the sediment samples. 

c<- read.table("euk_nosing_15s741r_EMBL3_nbotus.txt",h=T) #EMBL3 functional group by sample matrix with # OTUs
#columns 2:16 are the sediment samples. 

div <- read.csv("euk_nosing_15s741r_metrics_expanded.csv") #Exported metrics with lake info and approx dates added. 
#diversity metrics 

#Bind div to matrices a, b, and c and transform matrices. 
#OTUs - no. reads
a<- a[,1:16] #only need OTUs and samples
a<- melt(a) #make long-format
a<- dcast(a, variable ~ OTUs) #cast back into transposed matrix. 
a<- as.data.frame(cbind(div[,2:6], a[,2:884])) #bind with extra info
rownames(a)<- as.character(a[,2])

#Consider removing samples that don't meet a certain threshold for relative abundance 
a2<- a[,5:887]/rowSums(a[,5:887]) #matrix of OTU columns from a expressed as relative abundance for each sample. 

#0.01 - 1% relative abundance
summary(a2>0.01) 

a2.ind.01<- apply(a2, 2, function(x) any(x > 0.01)) #greater than 1% relative abundance 
a2.01.removals<- as.data.frame(a2[, !a2.ind.01]) #columns that don't contain a value greater than 0.01
a2.01.keep<- as.data.frame(a2[, a2.ind.01]) #columns that contain a value greater than 0.01
dim(a2.01.keep) #only 88 OTU columns 

#0.001 - 0.1% relative abundance 
summary(a2>0.001) 

a2.ind.001<- apply(a2, 2, function(x) any(x > 0.001)) #greater than 0.1% relative abundance 
a2.001.removals<- as.data.frame(a2[, !a2.ind.001]) 
a2.001.keep<- as.data.frame(a2[, a2.ind.001]) 
dim(a2.001.keep) #reduces to 883 species coumns, not a significant reduction 

#OTUs - no. reads, only OTUs with greater than 1% relative abundance in one sample 
a2.01<- as.data.frame(cbind(div[,2:5], a2.01.keep))

#EMBL3 - no. reads
b<- b[,1:16] #No. reads
b<- melt(b) 
b<- dcast(b, variable ~ EMBL3) 
b<- as.data.frame(cbind(div[,2:6], b[,2:37])) #bind with extra info
rownames(b)<- as.character(b[,2])

#EMBL3 - no. OTUs
c<- c[,1:16] #No. OTUs
c<- melt(c) 
c<- dcast(c, variable ~ EMBL3) 
c<- as.data.frame(cbind(div[,2:6], c[,2:37]))
rownames(c)<- as.character(c[,2])


#Unconstrained ordination (PCA)  #Keep in mind, width of data 
#PCA - sample by OTU (No. reads)
#Should probably screen out correlated OTUs
a.pca<- rda(a[,6:888]) #Should Hellinger transform reads- treat similar to counts?
a.sp.sc<- data.frame(scores(a.pca, choices = 1:2, display="sp"))
a.sit.sc<- data.frame(scores(a.pca, choices = 1:2, display="sites"))

a.pca.plot<-ggplot()
a.pca.plot<- a.pca.plot + geom_vline(x=0,colour="grey50") 
a.pca.plot<- a.pca.plot + geom_hline(y=0,colour="grey50") 
a.pca.plot<- a.pca.plot + geom_text(data = a.sit.sc, aes(x = PC1, y = PC2, label = a$Approx_yearCIC, colour = a$Lake, position = "jitter"),size = 7, angle = 0,  vjust = 1)
a.pca.plot<- a.pca.plot + scale_colour_manual(values = c("DAR" = "darkblue", "KB" = "red2"), "Lake")
a.pca.plot<- a.pca.plot + geom_point(data = a.sp.sc, aes(x = PC1, y = PC2), pch = 3, size = 2) 
a.pca.plot<- a.pca.plot + theme_bw()
a.pca.plot<- a.pca.plot + labs(x= "PC1 = <0.1% var exp.", y= "PC2 = <0.1%")
a.pca.plot<- a.pca.plot + theme(axis.text.x = element_text(colour="black", size=16))
a.pca.plot<- a.pca.plot + theme(axis.text.y = element_text(colour="black", size=16))
a.pca.plot<- a.pca.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
a.pca.plot<- a.pca.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
a.pca.plot<- a.pca.plot + ggtitle("PCA - OTUs (Reads non-transformed)")
#Add if want vectors
a.pca.plot<- a.pca.plot + geom_path(data = a.sit.sc, aes(x = PC1, y= PC2, colour= a$Lake), size =0.2, arrow=arrow(length= unit(0.4,"cm"), ends="first", type="open" ))


#PCA - sample by OTU (presence-absence)
a.pca<- rda(decostand(a[,6:888], method = "pa")) 
a.sp.sc<- data.frame(scores(a.pca, choices = 1:2, display="sp"))
a.sit.sc<- data.frame(scores(a.pca, choices = 1:2, display="sites"))

a.pca.plot<-ggplot()
a.pca.plot<- a.pca.plot + geom_vline(x=0,colour="grey50") 
a.pca.plot<- a.pca.plot + geom_hline(y=0,colour="grey50") 
a.pca.plot<- a.pca.plot + geom_text(data = a.sit.sc, aes(x = PC1, y = PC2, label = a$Approx_yearCIC, colour = a$Lake, position = "jitter"),size = 7, angle = 0,  vjust = 1)
a.pca.plot<- a.pca.plot + scale_colour_manual(values = c("DAR" = "darkblue", "KB" = "red2"), "Lake")
a.pca.plot<- a.pca.plot + geom_point(data = a.sp.sc, aes(x = PC1, y = PC2), pch = 3, size = 2) 
a.pca.plot<- a.pca.plot + theme_bw()
a.pca.plot<- a.pca.plot + labs(x= "PC1 = 19% var exp.", y= "PC2 = 17%")
a.pca.plot<- a.pca.plot + theme(axis.text.x = element_text(colour="black", size=16))
a.pca.plot<- a.pca.plot + theme(axis.text.y = element_text(colour="black", size=16))
a.pca.plot<- a.pca.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
a.pca.plot<- a.pca.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
a.pca.plot<- a.pca.plot + ggtitle("PCA - OTUs (Presence-Absence)")
#Add if want vectors
a.pca.plot<- a.pca.plot + geom_path(data = a.sit.sc, aes(x = PC1, y= PC2, colour= a$Lake), size =0.2, arrow=arrow(length= unit(0.4,"cm"), ends="first", type="open" ))


#PCA - sample by EMBL3 functional group (No. reads)
b.pca<- rda(b[,6:41]) 
b.sp.sc<- data.frame(scores(b.pca, choices = 1:2, display="sp"))
b.sit.sc<- data.frame(scores(b.pca, choices = 1:2, display="sites"))

b.pca.plot<-ggplot()
b.pca.plot<- b.pca.plot + geom_vline(x=0,colour="grey50") 
b.pca.plot<- b.pca.plot + geom_hline(y=0,colour="grey50") 
b.pca.plot<- b.pca.plot + geom_text(data = b.sit.sc, aes(x = PC1, y = PC2, label = b$Approx_yearCIC, colour = b$Lake, position = "jitter"),size = 7, angle = 0,  vjust = 1)
b.pca.plot<- b.pca.plot + scale_colour_manual(values = c("DAR" = "darkblue", "KB" = "red2"), "Lake")
b.pca.plot<- b.pca.plot + geom_text(data = b.sp.sc, aes(x = PC1, y = PC2, label = rownames(b.sp.sc)), size = 4) 
b.pca.plot<- b.pca.plot + theme_bw()
b.pca.plot<- b.pca.plot + labs(x= "PC1 = 79% var exp.", y= "PC2 = <1%")
b.pca.plot<- b.pca.plot + theme(axis.text.x = element_text(colour="black", size=16))
b.pca.plot<- b.pca.plot + theme(axis.text.y = element_text(colour="black", size=16))
b.pca.plot<- b.pca.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
b.pca.plot<- b.pca.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
b.pca.plot<- b.pca.plot + ggtitle("PCA - EMBL3 (Reads non-transformed)")
#Add if want vectors
b.pca.plot<- b.pca.plot + geom_path(data = b.sit.sc, aes(x = PC1, y= PC2, colour= a$Lake), size =0.2, arrow=arrow(length= unit(0.4,"cm"), ends="first", type="open" ))


#NMDS - sample by OTU (No. reads)
a.nmds<-metaMDS(a[,6:888],distance="bray") #warning messages associated with 'jaccard'
#Do with both jaccard and bray (Bray-Curtis) #same
a1 <- a.nmds$points[,1]
a2 <- a.nmds$points[,2]
ac<- as.data.frame(cbind(a1, a2))

a.nmds.plot<- ggplot()
a.nmds.plot<- a.nmds.plot + geom_vline(x=0,colour="grey50")
a.nmds.plot<- a.nmds.plot + geom_hline(y=0,colour="grey50")
a.nmds.plot<- a.nmds.plot + geom_text(data = ac, aes(x=a1, y=a2, label = a$Approx_yearCIC, colour = a$Lake), size = 7)
a.nmds.plot<- a.nmds.plot + scale_colour_manual(values = c("DAR" = "darkblue", "KB" = "red2"), "Lake")
a.nmds.plot<- a.nmds.plot + theme_bw()
a.nmds.plot<- a.nmds.plot + labs(x= "NMDS 1", y= "NMDS 2")
a.nmds.plot<- a.nmds.plot + theme(axis.text.x = element_text(colour="black", size=16))
a.nmds.plot<- a.nmds.plot + theme(axis.text.y = element_text(colour="black", size=16))
a.nmds.plot<- a.nmds.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
a.nmds.plot<- a.nmds.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
a.nmds.plot<- a.nmds.plot + ggtitle("NMDS - OTUs (Reads non-transformed)")
#Add if want vectors
a.nmds.plot<- a.nmds.plot + geom_path(data = ac, aes(x = a1, y= a2, colour= a$Lake), size =0.4, arrow=arrow(length= unit(0.4,"cm"), ends="first", type="open" ))


#NMDS - sample by OTU (No. reads) - only OTUs with at least 1% relative abundance in at least one sample
a.nmds<-metaMDS(a2.01[,5:92],distance="bray") 
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
a.nmds.plot<- a.nmds.plot + ggtitle("NMDS - OTUs (1%) (Reads non-transformed)")
#Add if want vectors
a.nmds.plot<- a.nmds.plot + geom_path(data = ac, aes(x = a1, y= a2, colour= a$Lake), size =0.4, arrow=arrow(length= unit(0.4,"cm"), ends="first", type="open" ))


#NMDS - sample by EMBL3 functional group (No. reads)
b.nmds<-metaMDS(b[,6:41],distance="bray") #warning messages associated with 'jaccard'
#Do with jaccard and bray. 
b1 <- b.nmds$points[,1]
b2 <- b.nmds$points[,2]
bc<- as.data.frame(cbind(b1, b2))

b.nmds.plot<- ggplot()
b.nmds.plot<- b.nmds.plot + geom_vline(x=0,colour="grey50")
b.nmds.plot<- b.nmds.plot + geom_hline(y=0,colour="grey50")
b.nmds.plot<- b.nmds.plot + geom_text(data = bc, aes(x=b1, y=b2, label = b$Approx_yearCIC, colour = b$Lake), size = 7)
b.nmds.plot<- b.nmds.plot + scale_colour_manual(values = c("DAR" = "darkblue", "KB" = "red2"), "Lake")
b.nmds.plot<- b.nmds.plot + theme_bw()
b.nmds.plot<- b.nmds.plot + labs(x= "NMDS 1", y= "NMDS 2")
b.nmds.plot<- b.nmds.plot + theme(axis.text.x = element_text(colour="black", size=16))
b.nmds.plot<- b.nmds.plot + theme(axis.text.y = element_text(colour="black", size=16))
b.nmds.plot<- b.nmds.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
b.nmds.plot<- b.nmds.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
b.nmds.plot<- b.nmds.plot + ggtitle("NMDS - EMBL3 (Reads non-transformed)")
#Add if want vectors
b.nmds.plot<- b.nmds.plot + geom_path(data = bc, aes(x = b1, y= b2, colour= b$Lake), size =0.4, arrow=arrow(length= unit(0.4,"cm"), ends="first", type="open" ))


#Hierarchical clustering tree - sample by OTU (No. reads)
a.bray<-vegdist(a[,6:888], method="bray")

a.clust<-hclust(a.bray) #surface of both cores cluster together.
plot(a.clust, hang=-1)

a.clust2<-reorder.hclust(a.clust, a.bray)
plot(a.clust2, hang=-1)

#Using ggdendro
ggdendrogram(a.clust2, rotate = FALSE, size = 2)

a.clust3 <- as.dendrogram(a.clust2)

clust.dat <- dendro_data(a.clust3, type = "rectangle")

clust.plot <- ggplot() + geom_segment(data = clust.dat$segments, aes(x = x, y = y, xend = xend, yend = yend)) 
clust.plot<- clust.plot + theme_bw()
clust.plot<- clust.plot + geom_text(data = clust.dat$labels, aes(x = x, y = y, label = a$Approx_yearCIC, colour = a$Lake), size = 3, vjust = 2)
clust.plot<- clust.plot + scale_colour_manual(values = c("DAR" = "darkblue", "KB" = "red2"), "Lake")


#Hierarchical clustering tree - sample by OTU (1% relative abundance)
a.bray<-vegdist(a2.01[,5:92], method="bray")

a.clust<-hclust(a.bray) #surface of both cores cluster together.
plot(a.clust, hang=-1)

a.clust2<-reorder.hclust(a.clust, a.bray)
plot(a.clust2, hang=-1)

#Using ggdendro
ggdendrogram(a.clust2, rotate = FALSE, size = 2)

a.clust3 <- as.dendrogram(a.clust2)

clust.dat <- dendro_data(a.clust3, type = "rectangle")

clust.plot <- ggplot() + geom_segment(data = clust.dat$segments, aes(x = x, y = y, xend = xend, yend = yend)) 
clust.plot<- clust.plot + theme_bw()
clust.plot<- clust.plot + geom_text(data = clust.dat$labels, aes(x = x, y = y, label = a$Approx_yearCIC, colour = a$Lake), size = 3, vjust = 2)
clust.plot<- clust.plot + scale_colour_manual(values = c("DAR" = "darkblue", "KB" = "red2"), "Lake")


##Superimpose clustering objects onto NMDS

#Sample by OTUs (No. reads) - 'a'
#Use NMDS and clustering objects run with matrix 'a' above, or run again. 
a.bray<-vegdist(a[,6:888], method="bray")
a.clust<-hclust(a.bray)
a.grp<- cutree(a.clust, 2) #cutting into 3 groups (?)- verify # of groups. 

a.nmds<-metaMDS(a[,6:888],distance="bray") 

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
a.bray<-vegdist(a[,6:888], method="bray") #distance matrix 
a.clust<-hclust(a.bray) #cluster
a.grp<- as.data.frame(cutree(a.clust, 2)) #define groups 
a.nmds<-metaMDS(a[,6:888],distance="bray") 
a1 <- a.nmds$points[,1]
a2 <- a.nmds$points[,2]
ac<- as.data.frame(cbind(a1, a2, a.grp))
colnames(ac)[3]<- 'grp'
grp<- as.factor(ac$grp)

#For manuscript figure (Fig. 4)
a.nmds.plot<- ggplot()
a.nmds.plot<- a.nmds.plot + geom_vline(x=0,colour="grey50")
a.nmds.plot<- a.nmds.plot + geom_hline(y=0,colour="grey50")
a.nmds.plot<- a.nmds.plot + geom_text(data = ac, aes(x=a1, y=a2, label = a$Approx_yearCIC, colour = grp), size = 7)
#a.nmds.plot<- a.nmds.plot + scale_colour_manual(values = c("1" = "darkblue", "2" = "red2"), "Cluster group")
a.nmds.plot<- a.nmds.plot + theme_bw()
a.nmds.plot<- a.nmds.plot + labs(x= "NMDS 1", y= "NMDS 2")
a.nmds.plot<- a.nmds.plot + theme(axis.text.x = element_text(colour="black", size=16))
a.nmds.plot<- a.nmds.plot + theme(axis.text.y = element_text(colour="black", size=16))
a.nmds.plot<- a.nmds.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
a.nmds.plot<- a.nmds.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
#a.nmds.plot<- a.nmds.plot + ggtitle("NMDS - OTUs (Reads non-transformed, with clusters)")
#Add if want vectors
a.nmds.plot<- a.nmds.plot + geom_path(data = ac, aes(x = a1, y= a2, line = a$Lake), size =0.4, arrow=arrow(length= unit(0.4,"cm"), ends="first", type="open" ))
a.nmds.plot<- a.nmds.plot + annotate("text", x=0.5, y=-0.5, label="Dauriat", size=10)
a.nmds.plot<- a.nmds.plot + annotate("text", x=-0.4, y=0, label="Knob", size=10)


#Sample by OTUs (1% relative abund, No. reads) - 'a2.01' 
#Use NMDS and clustering objects run with matrix 'a2.01' above, or run again. 
#Can come back to - NMDS wasn't too different when using reduced set of OTUs. 


#################MULTIVARIATE REGRESSION TREES##################
##Find breaks in OTU matrix in terms of explanatory variables

#Response matrices - should be read in, transformed and bound to 
#extra lake information the same way as for "Ordination analyses"
a #sample by OTUs (# of reads) [,6:888]
b #sample by EMBL3 (# of reads)[,6:41]
c #sample bu EMBL3 (# of OTUs) [,6:41]

#Explanatory variables - concern about too many for n = (check)
a$Approx_yearCIC  #year (same for b and c)

a$Time_period <- c("Post-mining", "Post-mining", "Post-mining", "Mining", "Mining", "Mining", "Pre-mining", "Pre-mining", "Pre-mining", "Post-mining", "Mining", "Post-mining", "Post-mining", "Post-mining", "Post-mining")
#Note assigned 1982 for KB to "mining"
a$Time_period2<- c(3,3,3,2,2,2,1,1,1,3,2,1,1,1,1) #dummy variable
#3 = "Post-mining", 2 = "Mining", 1 = "Pre-mining"
b$Time_period <- c("Post-mining", "Post-mining", "Post-mining", "Mining", "Mining", "Mining", "Pre-mining", "Pre-mining", "Pre-mining", "Post-mining", "Mining", "Post-mining", "Post-mining", "Post-mining", "Post-mining")
#Note assigned 1982 for KB to "mining"
b$Time_period2<- c(3,3,3,2,2,2,1,1,1,3,2,1,1,1,1) #dummy variable
#3 = "Post-mining", 2 = "Mining", 1 = "Pre-mining"
c$Time_period <- c("Post-mining", "Post-mining", "Post-mining", "Mining", "Mining", "Mining", "Pre-mining", "Pre-mining", "Pre-mining", "Post-mining", "Mining", "Post-mining", "Post-mining", "Post-mining", "Post-mining")
#Note assigned 1982 for KB to "mining"
c$Time_period2<- c(3,3,3,2,2,2,1,1,1,3,2,1,1,1,1) #dummy variable
#3 = "Post-mining", 2 = "Mining", 1 = "Pre-mining"

a$Lake #lake (same for b and c)- dummy code into 1 or 2
a$Lake2<- c(1,1,1,1,1,1,1,1,1,2,2,2,2,2,2) #DAR = 1, KB = 2
b$Lake2<- c(1,1,1,1,1,1,1,1,1,2,2,2,2,2,2) #DAR = 1, KB = 2
c$Lake2<- c(1,1,1,1,1,1,1,1,1,2,2,2,2,2,2) #DAR = 1, KB = 2

#Metal_EFs - extrapolated from geochemical data (from schefferville_lmerdata_MAR2015.xlsx)
met<- read.table("euk_metalEF.csv",h=T, sep=",")
a$met<- met$Metal_EF
b$met<- met$Metal_EF
c$met<- met$Metal_EF
  
#Use: $Approx_year2, $Time_period2, $Lake2, a$met



#Multivariate regression trees 

#mvpart - sample by OTUs (# of reads) (matrix a) - large model
y <- as.matrix(a[,5:887])
a.mrt<- mvpart(y ~ a$Approx_year2 + a$Time_period2 + a$Lake2 + a$met, data = a, cp=0, xv="pick", margin=0.08, xvmult=100, which=4) 
a.mrt<- mvpart(y ~ a$Approx_year2, data = a, cp=0, xv="pick", margin=0.08, xvmult=100, which=4)
#Need to ask about all these input paraemters
a.mrt
a.mrt$frame$dev
printcp(a.mrt) #Approx year and Lake retained. 
plotcp(a.mrt,minline=TRUE,upper=c("size"))
summary(a.mrt)
par(mfrow=c(1,2)) # two plots on one page 
rsq.rpart(a.mrt)  # visualize cross-validation result
#R2 is 1-Rel Error
# R2  
a.mrt$cptable
1-a.mrt$cptable[3,5]  

#Plot (MVPARTwrap pkg)
par(mfrow=c(1,1))
tree.a.mrt <- MRT(a.mrt,percent=10,species=colnames(y))
plot(tree.a.mrt)
summary(tree.a.mrt)
#Need to verify interpretation
#plot(as.party(tree.a.mrt), tp_args = list(id = FALSE)) #Not for MRT


#mvpart - sample by OTUs (# of reads - matrix a) - subset into individual lakes and
#just test with $Approx_year2 
a.DAR<- as.data.frame(subset(a, Lake == 'DAR', drop=T))
a.KB<- as.data.frame(subset(a, Lake == 'KB', drop=T))

y<- as.matrix(a.DAR[,6:888]) #swap out a.DAR with a.KB here and in mvpart() (next line)

a.mrt<- mvpart(y ~ a.DAR$Approx_yearCIC, data = a.DAR, cp=0, xv="pick", margin=0.08, xvmult=100, which=4) 
a.mrt
printcp(a.mrt) #Approx_yearCIC used
plotcp(a.mrt,minline=TRUE,upper=c("size"))
summary(a.mrt)
rsq.rpart(a.mrt)
a.mrt$cptable
1-a.mrt$cptable[2,3] #0.25 (DAR- use 3,3) #KB use (2,3) = 0.89
#Plot
tree.a.mrt <- MRT(a.mrt,percent=10,species=colnames(y))
plot(tree.a.mrt)
summary(tree.a.mrt)

#Test for $met splits in tree with both lakes
y<- as.matrix(a[,5:887]) #both lakes

a.mrt<- mvpart(y ~ a$met, data = a, cp=0, xv="pick", margin=0.08, xvmult=100, which=4) 
a.mrt
printcp(a.mrt) #
plotcp(a.mrt,minline=TRUE,upper=c("size"))
summary(a.mrt)
rsq.rpart(a.mrt)
a.mrt$cptable
1-a.mrt$cptable[4,3] #
#Plot
tree.a.mrt <- MRT(a.mrt,percent=10,species=colnames(y))
plot(tree.a.mrt)
summary(tree.a.mrt)

#Make a plot of the metal EFs in order to visualize where the breaks are
met.plot<- ggplot(a, aes(x=Approx_year2, y=met, fill=Lake)) + geom_bar(stat = 'identity', width = 5.00) + coord_flip() + facet_wrap(~Lake)
met.plot<- met.plot + scale_fill_manual(values = c("DAR" = "darkblue", "KB" = "red2"))
met.plot<- met.plot + labs(x= "Year", y= "Metal EF") + theme_bw()
met.plot<- met.plot + theme(axis.text.x = element_text(colour="black", size=16))
met.plot<- met.plot + theme(axis.text.y = element_text(colour="black", size=16))
met.plot<- met.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
met.plot<- met.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
#could annotate with segments, but not sure how to annotate in separate facets?
#http://stackoverflow.com/questions/11889625/annotating-text-on-individual-facet-in-ggplot2



#mvpart - sample by EMBL3 (# of reads) (matrix b)
y <- as.matrix(b[,5:40])
b.mrt<- mvpart(y ~ b$Approx_year2 + b$Time_period2 + b$Lake2 + b$met, data = b, cp=0, xv="pick", margin=0.08, xvmult=100, which=4) 
#Need to ask about all these input paraemters
b.mrt
b.mrt$frame$dev
printcp(b.mrt) #Approx year and Lake retained. 
plotcp(b.mrt,minline=TRUE,upper=c("size"))
summary(b.mrt)
#par(mfrow=c(1,2)) # two plots on one page 
rsq.rpart(b.mrt)  # visualize cross-validation result
#R2 is 1-Rel Error
# R2  
b.mrt$cptable
1-b.mrt$cptable[5,3]  

#Plot (MVPARTwrap pkg)
tree.b.mrt <- MRT(b.mrt,percent=10,species=colnames(y))
plot(tree.b.mrt)
summary(tree.b.mrt)


#mvpart - sample by EMBL3 (# of OTUs) (matrix c)
y <- as.matrix(c[,5:40])
c.mrt<- mvpart(y ~ c$Approx_year2 + c$Time_period2 + c$Lake2 + c$met, data = c, cp=0, xv="pick", margin=0.08, xvmult=100, which=4) 
#Need to ask about all these input paraemters
c.mrt
c.mrt$frame$dev
printcp(c.mrt) #Approx year and Lake retained. 
plotcp(c.mrt,minline=TRUE,upper=c("size"))
summary(c.mrt)
#par(mfrow=c(1,2)) # two plots on one page 
rsq.rpart(c.mrt)  # visualize cross-validation result
#R2 is 1-Rel Error
# R2  
c.mrt$cptable
1-c.mrt$cptable[3,5]  

#Plot (MVPARTwrap pkg)
tree.c.mrt <- MRT(c.mrt,percent=10,species=colnames(y))
plot(tree.c.mrt)
summary(tree.c.mrt)


#################LINEAR MIXED EFFECT MODELS##################
##Test the relationship between Eukaryote Shannon diversity and 
#Metal_EF, consideirng lake and time period as random factors. 

div<- read.csv("euk_nosing_15s741r_metrics_expanded.csv")
met<- read.table("euk_metalEF.csv",h=T, sep=",") 

d<- as.data.frame(cbind(div, met$Metal_EF))
colnames(d) [16]<- 'Metal_EF'
d$Time_period <- c("Post-mining", "Post-mining", "Post-mining", "Mining", "Mining", "Mining", "Pre-mining", "Pre-mining", "Pre-mining", "Post-mining", "Mining", "Post-mining", "Post-mining", "Post-mining", "Post-mining")


#Null- H0 - linear model 
lm1<- lm(Shannon~Metal_EF, data=d) 
summary(lm1) #Adj R2 = 0.123
lm1.resid<- rstandard(lm1)
#Lake effect
plot(lm1.resid~ d$Lake, xlab = "Lake", ylab="Standardized residuals")
abline(0,0, lty=2) #show variation across lakes 
#Time period effect
plot(lm1.resid~ d$Time_period, xlab = "Time_period", ylab="Standardized residuals")
abline(0,0, lty=2) #issue

##Mixed models with random factors are 
#appropriate for all the models- are they appropriate for time period? 

#Base Model 1: Euk_Shannon~Metal_EF + Lake + Time_period
#Note REML = TRUE for these models 

#Base 1 model 1
#Description: Lake and time period are random factors, varying intercepts for both
b1m1<- lmer(Shannon ~ Metal_EF + (1|Lake) + (1|Time_period), data=d, REML=TRUE)
b1m1
summary(b1m1)

#Base 1 model 2
#Description: Lake and time period as random factors, intercepts and slopes vary with respect to metal
b1m2<- lmer(Shannon~Metal_EF + (1+Metal_EF|Lake) + (1+Metal_EF|Time_period), data=d, REML=TRUE)
b1m2
summary(b1m2)

#Base 1 model 3
#Description: Lake as a random factor, varying intercept 
b1m3<- lmer(Shannon~Metal_EF + (1|Lake), data=d, REML=TRUE)
b1m3
summary(b1m3)

#Base 1 model 4
#Description: Time period as random factor, varying intercept 
b1m4<- lmer(Shannon~Metal_EF + (1|Time_period), data=d, REML=TRUE)
b1m4
summary(b1m4)

#Base 1 model 5
#Description: Lake as a random factor, varying slope and intercept with respect to metal 
b1m5<- lmer(Shannon~Metal_EF + (1+Metal_EF|Lake), data=d, REML=TRUE)
b1m5
summary(b1m5)

#Base 1 model 6
#Description: Time period as a random factor, varying slope and intercept with respect to metal 
b1m6<- lmer(Shannon~Metal_EF + (1+Metal_EF|Time_period), data=d, REML=TRUE)
b1m6
summary(b1m6)

#Base 1 model 7
#Description: Lake and time period as random factors, varying intercept plus slope for lake, varying intercept only for time period 
b1m7<- lmer(Shannon~Metal_EF + (1+Metal_EF|Lake) + (1|Time_period), data=d, REML=TRUE)
b1m7
summary(b1m7)

#Base 1 model 8
#Description: Lake and time period as random factors, varying intercept plus slope for time period, varying intercept only for lake
b1m8<- lmer(Shannon~Metal_EF + (1|Lake) + (1+Metal_EF|Time_period), data=d, REML=TRUE)
b1m8
summary(b1m8) #Model fail

##Model selection 
#Note REML = FALSE for model selection/comparison  

b1m1b<- lmer(Shannon ~ Metal_EF + (1|Lake) + (1|Time_period), data=d, REML=FALSE)

b1m2b<- lmer(Shannon~Metal_EF + (1+Metal_EF|Lake) + (1+Metal_EF|Time_period), data=d, REML=FALSE)

b1m3b<- lmer(Shannon~Metal_EF + (1|Lake), data=d, REML=FALSE)

b1m4b<- lmer(Shannon~Metal_EF + (1|Time_period), data=d, REML=FALSE)

b1m5b<- lmer(Shannon~Metal_EF + (1+Metal_EF|Lake), data=d, REML=FALSE)

b1m6b<- lmer(Shannon~Metal_EF + (1+Metal_EF|Time_period), data=d, REML=FALSE)

b1m7b<- lmer(Shannon~Metal_EF + (1+Metal_EF|Lake) + (1|Time_period), data=d, REML=FALSE)

lm1<- lm(Shannon~Metal_EF, data=d)

##Model evaluation
AICc<-c(AICc(b1m1b), AICc(b1m2b), AICc(b1m3b), AICc(b1m4b), AICc(b1m5b), AICc(b1m6b), AICc(b1m7b), AICc(lm1))
# Put values into one table for easy comparision
Model<-c("b1m1b", "b1m2b", "b1m3b", "b1m4b", "b1m5b", "b1m6b", "b1m7b", "lm1")
AICtable<-data.frame(Model=Model, AICc=AICc)
AICtable #lowest AIC is the null linear model. 

##Continue investigation (scheff_LMER_March2015.R script)


#################ANOSIM##################
##Statistically test whether there is a signifcant difference
#between groups (vegan pkg)

#Using anosim- investigate adonis.****
anosim(dat, grouping, permutations = 999, distance = "bray", strata)

#Test significant differences between lakes - OTUs (reads)
#matrix a from MRT analyses. 
a.ano<- anosim(a[,5:887], a$Lake, permutations = 999, distance = "bray")
summary(a.ano)

#Test significant differences between time periods - OTUs (reads)
#matrix a from MRT analyses. 
a.ano<- anosim(a[,5:887], a$Time_period, permutations = 999, distance = "bray")
summary(a.ano)

#Test time period effect for lakes separately 
a.DAR<- as.data.frame(subset(a, Lake == "DAR", drop=T))
a.KB<- as.data.frame(subset(a, Lake == "KB", drop=T))

a.DAR.ano<- anosim(a.DAR[,5:887], a.DAR$Time_period, permutations = 999, distance = "bray")
summary(a.DAR.ano)

a.KB.ano<- anosim(a.KB[,5:887], a.KB$Time_period, permutations = 999, distance = "bray")
summary(a.KB.ano)

#Test if the 2 groups identified in the hierarchical clustering are sig different
a$clustgrp<- c(1,1,1,1,2,2,1,2,2,1,1,2,2,2,2)
as.factor(a$clustgrp)

a.ano<- anosim(a[,5:887], a$clustgrp, permutations = 999, distance = "bray")
summary(a.ano)

#Test significant differences between lakes - EMBL3 (reads)
#matrix b from MRT analyses. 
b.ano<- anosim(b[,5:40], b$Lake, permutations = 999, distance = "bray")
summary(b.ano)

#Test significant differences between time periods - EMBL3 (reads)
#matrix b from MRT analyses. 
b.ano<- anosim(b[,5:40], b$Time_period, permutations = 999, distance = "bray")
summary(b.ano)

#Test significant differences between lakes - EMBL3 (# OTUs)
#matrix c from MRT analyses. 
c.ano<- anosim(c[,5:40], c$Lake, permutations = 999, distance = "bray")
summary(c.ano)

#Test significant differences between time periods - EMBL3 (# OTUs)
#matrix c from MRT analyses. 
c.ano<- anosim(c[,5:40], c$Time_period, permutations = 999, distance = "bray")
summary(c.ano)


#################TEMPORAL BETA-DIVERSITY##################
##Use Legendre functions to compute temporal beta diversity
#between time periods and years. 

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

##paired.diff2 function 
#####
paired.diff2 <- function(mat1,mat2,method="hellinger", pa.tr=FALSE, nperm=99, permute.sp=1, BCD=TRUE, replace=FALSE, clock=FALSE)
  #
  # Compute and test differences between pairs of data vectors of observations at T1 and T2.
  #
  # Test hypothesis (H0) that an object is not exceptionally different between T1 and T2.
  # Example in palaeoecology: ancient and modern diatom communities in sediment cores.
  # Example in sequence data: "regions" before and after treatment.
  # 
  # Arguments --
  # mat1, mat2: two matrices (or data.frames) with the same number of rows and columns. 
  # The rows must correspond to the same objects (e.g. sites) and the colums to the same 
  # variables (e.g. species).
# method={"hellinger", "chord", "ruzicka", "%difference", "euclidean"}. 
#   Methods {"hellinger", "chord"} are obtained by transformation of the  
#     species data followed by calculation of the Euclidean distance. These distances 
#     have the Euclidean property. 
#     If pa.tr=TRUE, sqrt(2)*sqrt(1-Ochiai) is computed.
#   Methods {"ruzicka", "%difference"} are obtained by computing a dissimilarity function. 
#     It is recommended to take the square root of these dissimilarities before computing 
#     ordinations by principal coordinate analysis. However, that precaution is not 
#     important here; the results of the permutation tests will be the same for these
#     dissimilarities square-rooted or not.
#     If pa.tr=TRUE, either the Jaccard or the Sørensen coefficient is computed.
# pa.tr=FALSE: do NOT transform the data to presence-absence.
#      =TRUE : transform the data to binary (i.e. presence-absence) form.
# nperm = number of permutations for the permutation test.
##
# This version of the function contains three permutation methods --
# permute.sp=1 : permute data separately in each column, both matrices in the same way.
#           =2 : permute data separately in each column. Do not force the permutations to 
#  			 start at the same point in the two matrices.
#           =3 : permute entire rows in each matrix separately (suggestion D. Borcard).
##
# BCD=TRUE  : Compute and save the B and C components of the %difference and Ruzicka D.
#             For the %difference, they are expressed as B/(2A+B+C) and C/(2A+B+C).
#             For the Ruzicka D, they are expressed as B/(A+B+C) and C/(A+B+C).
# BCD=FALSE : Do not compute the components. BCD=FALSE for D other than %diff and Ruzicka.
# replace=FALSE : sampling without replacement for regular permutation test.
#        =TRUE  : sampling with replacement. The testing method is then bootstrapping.
# clock=FALSE : Do not print the computation time.
#      =TRUE  : Print the time (in sec) used for computation.
#
# Details --
# H0: in each matrix (e.g. each time), the sites do not differ in species composition. 
#     They only differ by random sampling of each species' statistical population.
# H1: Some sites are exceptionally different between T1 and T2.
# 
# The randomization procedures are the following:
# 1. In each matrix, the original values (e.g. species abundances) are permuted at random, 
# independently in each column. Permutation of the two matrices is started with the same 
# random seed, so that the values in each column (e.g. species) are permuted in the same 
# way in mat1.perm and mat2.perm. 
# 2. The transformation, if any, is recomputed on the permuted data matrices. This is 
# necessary to make sure that the permuted data are transformed in the same way as the 
# initial data, with row sums or row lengths of 1. In this way, the D of the permuted data 
# will be comparable to the reference D.
# 3. The distances between T1 and T2 are recomputed, for each site separately.
#
# For presence-absence data, this function computes the binary forms of the quantitative 
# coefficients listed under the 'method' parameter. The "hellinger" and "chord" 
# transformations produce the Ochiai distance, or more precisely: 
# D.Hellinger = D.chord = sqrt(2) * sqrt(1 - S.Ochiai) 
# where "S.Ochiai" designates the Ochiai similarity coefficient.  
# The "%difference" dissimilarity produces (1 – S.Sørensen) 
# whereas the "ruzicka" dissimilarity produces (1 – S.Jaccard).
#
# Community composition data could be log-transformed prior to analysis. Only the 
# Euclidean distance option should be used with log-transformed data. It is meaningless to 
# subject log-transformed data to the {"hellinger", "chord"} transformations 
# available in this function. - One can use either the log(y+1 transformation (log1p() 
# function of {base}), or Anderson et al. (2006) special log transformation available in 
# {vegan}: decostand(mat, "log", logbase=10).
#
# Value --
# A list containing the vector of distances between T1 and T2 for each object and a 
# corresponding vector of p-values. The significant p-values (e.g. p.dist ≤ 0.05) indicate 
# exceptional objects for the difference of their species composition. The p-values should be corrected for multiple testing using function p.adjust() of {stats}. A good general choice is method="holm", which is the default option of the function.
# An output table containing B, C and D.
#
# Author:: Pierre Legendre
# License: GPL-2 
{
  ### Internal functions
  dissim <- function(mat1, mat2, n, method, tr=TRUE, BCD, ref)
  {
    vecD = vector(mode="numeric",length=n)
    if(BCD) { 
      vecB = vector(mode="numeric",length=n)
      vecC = vector(mode="numeric",length=n)
      vecD = vector(mode="numeric",length=n)
    } else { vecB=NA; vecC=NA; vecD=NA }
    #
    # Compute the dissimilarity between T1 and T2 for each object (site)
    # 1. If method = {"hellinger", "chord"}, tr is TRUE
    if(tr) for(i in 1:n) vecD[i] = dist(rbind(mat1[i,], mat2[i,]))
    #
    # 2. Compute the Euclidean distance
    if(method == "euclidean")  
      for(i in 1:n) vecD[i] = dist(rbind(mat1[i,], mat2[i,])) 
      # 3. Compute the Ruzicka or %difference dissimilarity 
      if(method == "ruzicka") dissimil=1       # Quantitative form of Jaccard
      if(method == "%difference") dissimil=2   # Quantitative form of Sørensen
      if(any(method == c("ruzicka", "%difference"))) { 
        for(i in 1:n) {
          tmp = RuzickaD(mat1[i,], mat2[i,], method=method, BCD=BCD, ref=ref) 
          if(BCD) {
            vecB[i] <- tmp$B
            vecC[i] <- tmp$C }
          vecD[i] <- tmp$D
        }
      }
      # Alternative method (not used here) to compute the %difference dissimilarity:
      #	for(i in 1:n) vecD[i] = vegdist(rbind(mat1[i,], mat2[i,]), "bray")         #Slower
      list(vecB=vecB, vecC=vecC, vecD=vecD)
  }
  ###
  transform <- function(mat, method)
  {
    if(method=="hellinger") mat = decostand(mat, "hellinger")
    if(method=="chord")     mat = decostand(mat, "norm")
    mat
  }
  ### End internal functions
  ###
  A <- system.time({
    
    epsilon <- sqrt(.Machine$double.eps)
    method <- match.arg(method, c("euclidean", "hellinger", "chord", "ruzicka", "%difference")) 
    n = nrow(mat1)
    p = ncol(mat1)
    if((nrow(mat2)!=n) | (ncol(mat2)!=p)) stop("The matrices are not of the same size.")
    #
    if(pa.tr) {
      mat1 <- ifelse(mat1>0, 1, 0)
      mat2 <- ifelse(mat2>0, 1, 0) }
    if(any(method == c("hellinger", "chord"))) {
      tr <- TRUE
      require(vegan)
    } else { tr <- FALSE }
    if( (any(method == c("ruzicka", "%difference"))) & BCD) { 
      BCD.mat <- matrix(0,n,3)
      if(method=="ruzicka")    colnames(BCD.mat) <- 
          c("B/(A+B+C)","C/(A+B+C)","D=(B+C)/(A+B+C)")
      if(method=="%difference") colnames(BCD.mat) <- 
          c("B/(2A+B+C)","C/(2A+B+C)","D=(B+C)/(2A+B+C)")
      rownames(BCD.mat) <- paste("Obj",1:n,sep=".")
    } else {
      BCD <- FALSE 
      BCD.mat <- NA }
    ###
    # 1. Compute the reference D for each object from corresponding vectors in the 2 matrices.
    if(tr) { 
      tmp <-dissim(transform(mat1,method), transform(mat2,method),n,method,tr,BCD,ref=FALSE)
    } else { tmp <- dissim(mat1, mat2, n, method, tr, BCD, ref=TRUE) }
    vecD.ref <- tmp$vecD
    if(BCD) { BCD.mat[,1]<-tmp$vecB ; BCD.mat[,2]<-tmp$vecC ; BCD.mat[,3]<-tmp$vecD }
    ###
    if(permute.sp!=3) {   # Permute the data separately in each column.
      # 2. Permutation methods 1 and 2 --
      # Permute *the raw data* by columns. Permute the two matrices in the same way, saving the seed before the two sets of permutations through sample(). 
      # Permutation test for each distance in vector D
      # seed: seed for random number generator, used by the permutation function 
      #       sample(). It is reset to that same value before permuting the values in the  
      #       columns of the second matrix. 
      if(nperm>0) {
        nGE.D = rep(1,n)
        for(iperm in 1:nperm) {
          BCD <- FALSE
          if(permute.sp==1) {    # Permutation methods 1
            seed <- ceiling(runif(1,max=100000))
            # cat("seed =",seed,'\n')
            set.seed(seed)
            mat1.perm <- apply(mat1,2,sample)
            set.seed(seed)
            mat2.perm <- apply(mat2,2,sample)
          } else {  # Permutation methods 2 - Do not force the permutations 
            # to start at the same point in the two matrices.
            mat1.perm <- apply(mat1,2,sample)
            mat2.perm <- apply(mat2,2,sample)
          }
          # 3. Recompute transformations of the matrices and the D values of the paired vectors.
          if(tr) { tmp <- dissim(transform(mat1.perm,method), 
                                 transform(mat2.perm,method), n, method, tr, BCD, ref=FALSE)
          } else { tmp <- dissim(mat1.perm, mat2.perm, n, method, tr, BCD, ref=FALSE) }
          vecD.perm <- tmp$vecD
          ge <- which(vecD.perm+epsilon >= vecD.ref)
          nGE.D[ge] <- nGE.D[ge] + 1
        }
        # 4. Compute the p-value associated with each distance.
        p.dist <- nGE.D/(nperm+1)
      } else { p.dist <- NA }   # if nperm=0
      
    } else if(permute.sp==3) {   
      # 2.bis  Permutation method 3 -- 
      # Permute entire rows in each matrix separately.
      if(nperm>0) {
        seed <- ceiling(runif(1,max=100000))
        set.seed(seed)
        nGE.D = rep(1,n)
        for(iperm in 1:nperm) {
          BCD <- FALSE
          mat1.perm <- mat1[sample(n),]
          mat2.perm <- mat2[sample(n),]
          #
          # 3.bis Recompute the D values of the paired vectors.
          if(tr) { tmp <- dissim(transform(mat1.perm,method), 
                                 transform(mat2.perm,method), n, method, tr, BCD, ref=FALSE)
          } else { tmp <- dissim(mat1.perm, mat2.perm, n, method, tr, BCD, ref=FALSE) }
          vecD.perm <- tmp$vecD
          ge <- which(vecD.perm+epsilon >= vecD.ref)
          nGE.D[ge] <- nGE.D[ge] + 1
        }
        # 4.bis Compute the p-value associated with each distance.
        p.dist <- nGE.D/(nperm+1)
      } else { p.dist <- NA }   # if nperm=0
    }
    p.adj <- p.adjust(p.dist,"holm")
  })
  A[3] <- sprintf("%2f",A[3])
  if(clock) cat("Computation time =",A[3]," sec",'\n')
  #
  list(vecD.ref=vecD.ref, p.dist=p.dist, p.adj=p.adj, BCD.mat=BCD.mat)
}

RuzickaD <- function(vec1, vec2, method="ruzicka", BCD=FALSE, ref=TRUE)
  #
  # Compute the Ruzicka dissimilarity (quantitative form of the Jaccard dissimilarity)
  # or the percentage difference (quantitative form of the Sørensen dissimilarity).
  # A single dissimilarity is computed because there are only two vectors in this function.
  #
  # Arguments --
  # vec1, vec2 : data vectors (species abundance or presence-absence data)
  # method == c("ruzicka", "%difference")
  # BCD=TRUE  : Compute and save the B and C components of the %difference and Ruzicka D.
  #             For the %difference, they are B/(2A+B+C), C/(2A+B+C), D/(2A+B+C).
  #             For the Ruzicka D, they are B/(A+B+C), C/(A+B+C), D/(A+B+C).
# BCD=FALSE : Do not compute the components. BCD=FALSE for D other than %diff and Ruzicka.
# ref=TRUE  : Compute the reference values of D, B and C
#    =FALSE : Under permutation, compute only the value of D. Use separate code (shorter).
#
# License: GPL-2 
# Author:: Pierre Legendre, April 2015
{
  # An algorithm applicable to matrices Y containing two data vectors only
  #
  A <- sum(pmin(vec1, vec2))          # A = sum of minima from comparison of the 2 vectors
  sum.Y <- sum(vec1, vec2)            # Sum of all values in the two vectors, (2A+B+C)
  #
  if(ref) {    # Compute the reference values of statistics D, B and C
    tmp = vec1 - vec2
    B = sum(tmp[tmp>0])                 # Sum of the species losses between T1 and T2
    C = -sum(tmp[tmp<0])                # Sum of the species gains between T1 and T2
    D = B+C                             # Dissimilarity
    
    # Under permutation, compute only the value of D. - Shorter computation time.
  } else { 
    D <- sum.Y-2*A                      # (B+C)
  }
  # Compute the denominator (den) of the Ruzicka or %difference index
  if(method == "ruzicka") { den <-(sum.Y-A)  # den = (A+B+C)
  } else { den <- sum.Y }                # den = (2A+B+C)
  if(!BCD) { B <- NA ; C <- NA }
  list(B=B/den, C=C/den, D=D/den)
}

# Examples -- 
# data(mite)
# res1 = paired.diff(mite[1:10,],mite[61:70,],method="hellinger",nperm=999,permute.sp=1)
# Computation time = 1.971000  sec 
# res2 = paired.diff(mite[1:10,],mite[61:70,],method="hellinger",nperm=999,permute.sp=3)


#####Need to use RAW data (non-rarefied for these analyses)
##READ IN RAW DATA AND ADD DESCRIPTIVE COLUMNS ETC. *****ALL TEMPORAL BD WITH NON-RAREFIED DATA
a <- read.table("euk_nosing_ss_EMBL.csv",h=T, sep=",") #eukaryote data, no singletons, matched to taxonomy
#Remove column 11 (Sed10)
a<- a[,-11]

div <- read.csv("euk_nosing_15s741r_metrics_expanded.csv") #For descriptive columns

#Add descriptive columns
a<- a[,1:16] #only need OTUs and samples
a<- melt(a) #make long-format
a<- dcast(a, variable ~ OTUs) #cast back into transposed matrix. 
a<- as.data.frame(cbind(div[,2:6], a[,2:884])) #bind with extra info

#Subset into the two lakes
a.DAR<- as.data.frame(subset(a, Lake == "DAR", drop=T))
a.KB<- as.data.frame(subset(a, Lake == "KB", drop=T))


##DAR - OTUs (reads) 
#Make each row a dataframe
#T0 - DAR9 (oldest) - 1930
dar.intT0<- as.data.frame(subset(a.DAR, Sample_ID == "DAR9", drop=T))

#T1 - DAR8 - 1945
dar.intT1<- as.data.frame(subset(a.DAR, Sample_ID == "DAR8", drop=T))

#T2 - DAR7 - 1950
dar.intT2<- as.data.frame(subset(a.DAR, Sample_ID == "DAR7", drop=T))

#T3 - DAR6 - 1963
dar.intT3<- as.data.frame(subset(a.DAR, Sample_ID == "DAR6", drop=T))

#T4 - DAR5 - 1973
dar.intT4<- as.data.frame(subset(a.DAR, Sample_ID == "DAR5", drop=T))

#T5 - DAR4 - 1981
dar.intT5<- as.data.frame(subset(a.DAR, Sample_ID == "DAR4", drop=T))

#T6 - DAR3 - 1994
dar.intT6<- as.data.frame(subset(a.DAR, Sample_ID == "DAR3", drop=T))

#T7 - DAR2 - 2004
dar.intT7<- as.data.frame(subset(a.DAR, Sample_ID == "DAR2", drop=T))

#T8 - DAR1 - 2011
dar.intT8<- as.data.frame(subset(a.DAR, Sample_ID == "DAR1", drop=T))

#Comparisons
#T0-T1 - 1930-1945
bdtemp.darT0T1.int<- decompose.D2(dar.intT0[,6:888], dar.intT1[,6:888], den.type=2)  
bdtemp.darT0T1.int.mat2<- as.data.frame(bdtemp.darT0T1.int$mat2) #making result matrix into a dataframe

#T1-T2 - 1945-1950
bdtemp.darT1T2.int<- decompose.D2(dar.intT1[,6:888], dar.intT2[,6:888], den.type=2)  
bdtemp.darT1T2.int.mat2<- as.data.frame(bdtemp.darT1T2.int$mat2) 

#T2-T3 - 195-1963
bdtemp.darT2T3.int<- decompose.D2(dar.intT2[,6:888], dar.intT3[,6:888], den.type=2)  
bdtemp.darT2T3.int.mat2<- as.data.frame(bdtemp.darT2T3.int$mat2) 

#T3-T4 - 1963-1973
bdtemp.darT3T4.int<- decompose.D2(dar.intT3[,6:888], dar.intT4[,6:888], den.type=2)  
bdtemp.darT3T4.int.mat2<- as.data.frame(bdtemp.darT3T4.int$mat2) 

#T4-T5 - 1973-1981
bdtemp.darT4T5.int<- decompose.D2(dar.intT4[,6:888], dar.intT5[,6:888], den.type=2)  
bdtemp.darT4T5.int.mat2<- as.data.frame(bdtemp.darT4T5.int$mat2) 

#T5-T6 - 1981-1994
bdtemp.darT5T6.int<- decompose.D2(dar.intT5[,6:888], dar.intT6[,6:888], den.type=2)  
bdtemp.darT5T6.int.mat2<- as.data.frame(bdtemp.darT5T6.int$mat2) 

#T6-T7 - 1994-2004
bdtemp.darT6T7.int<- decompose.D2(dar.intT6[,6:888], dar.intT7[,6:888], den.type=2)  
bdtemp.darT6T7.int.mat2<- as.data.frame(bdtemp.darT6T7.int$mat2) 

#T7-T8 - 2004-2010
bdtemp.darT7T8.int<- decompose.D2(dar.intT7[,6:888], dar.intT8[,6:888], den.type=2)  
bdtemp.darT7T8.int.mat2<- as.data.frame(bdtemp.darT7T8.int$mat2) 


#Between time periods - NOT SURE THIS IS APPROPRIATE
#Subset DAR into the three "Time_period_grouping" 
#a.DAR.PostMin<- as.data.frame(subset(a.DAR, Time_period == "Post-mining", drop=T))
#a.DAR.Min<- as.data.frame(subset(a.DAR, Time_period == "Mining", drop=T))
#a.DAR.PreMin<- as.data.frame(subset(a.DAR, Time_period == "Pre-mining", drop=T))

#PostMin - take average within time period
#a.DAR.PostMin.avg<- as.data.frame(colMeans(a.DAR.PostMin[,5:887]))
#a.DAR.PostMin.avg<- as.data.frame(t(a.DAR.PostMin.avg))

#Min - take average within time period
#a.DAR.Min.avg<- as.data.frame(colMeans(a.DAR.PostMin[,5:887]))
#a.DAR.Min.avg<- as.data.frame(t(a.DAR.Min.avg))

#PreMin - take average within time period
#a.DAR.PreMin.avg<- as.data.frame(colMeans(a.DAR.PreMin[,5:887]))
#a.DAR.PreMin.avg<- as.data.frame(t(a.DAR.PreMin.avg))

#Temporal BD between the different time periods
#PreMin to Min
#bdtemp.darPre_Min<- decompose.D2(a.DAR.PreMin.avg, a.DAR.Min.avg, den.type=2)  
#Min to PostMin
#bdtemp.darMin_Post<- decompose.D2(a.DAR.Min.avg, a.DAR.PostMin.avg, den.type=2)  
#PreMin to PostMin
#bdtemp.darPre_Post<- decompose.D2(a.DAR.PreMin.avg, a.DAR.PostMin.avg, den.type=2)  

#Looking at output of total temporal beta diversity and species gain/loss components 
#Pre-mining to mining
#darPre_Min.mat2<- as.data.frame(bdtemp.darPre_Min$mat2)

#Mining to post-mining
#darMin_Post.mat2<- as.data.frame(bdtemp.darMin_Post$mat2)

#Pre-mining to post-mining
#darPre_Post.mat2<- as.data.frame(bdtemp.darPre_Post$mat2)


##KB
#Make each row a dataframe
#T0 - KB6 (oldest) - Pre-1850
kb.intT0<- as.data.frame(subset(a.KB, Sample_ID == "KB6", drop=T))

#T1 - KB5 - Pre-1850
kb.intT1<- as.data.frame(subset(a.KB, Sample_ID == "KB5", drop=T))

#T2 - KB4 - 1866
kb.intT2<- as.data.frame(subset(a.KB, Sample_ID == "KB4", drop=T))

#T3 - KB3 - 1928
kb.intT3<- as.data.frame(subset(a.KB, Sample_ID == "KB3", drop=T))

#T4 - KB2 - 1982
kb.intT4<- as.data.frame(subset(a.KB, Sample_ID == "KB2", drop=T))

#T5 - KB1 (2012)
kb.intT5<- as.data.frame(subset(a.KB, Sample_ID == "KB1", drop=T))

#Comparisons
#T0-T1- Pre1850(1) to Pre1850(2)
bdtemp.kbT0T1.int<- decompose.D2(kb.intT0[,5:887], kb.intT1[,5:887], den.type=2)  
bdtemp.kbT0T1.int.mat2<- as.data.frame(bdtemp.kbT0T1.int$mat2) 

#T1-T2- Pre1850(2) to 1866
bdtemp.kbT1T2.int<- decompose.D2(kb.intT1[,5:887], kb.intT2[,5:887], den.type=2)  
bdtemp.kbT1T2.int.mat2<- as.data.frame(bdtemp.kbT1T2.int$mat2) 

#T2-T3- 1866 ro 1928
bdtemp.kbT2T3.int<- decompose.D2(kb.intT2[,5:887], kb.intT3[,5:887], den.type=2)  
bdtemp.kbT2T3.int.mat2<- as.data.frame(bdtemp.kbT2T3.int$mat2) 

#T3-T4- 1928 to 1982
bdtemp.kbT3T4.int<- decompose.D2(kb.intT3[,5:887], kb.intT4[,5:887], den.type=2)  
bdtemp.kbT3T4.int.mat2<- as.data.frame(bdtemp.kbT3T4.int$mat2) 

#T4-T5- 1982 to 2012
bdtemp.kbT4T5.int<- decompose.D2(kb.intT4[,5:887], kb.intT5[,5:887], den.type=2)  
bdtemp.kbT4T5.int.mat2<- as.data.frame(bdtemp.kbT4T5.int$mat2) 

#Between time period #Aren't really able to do comparisons between
#time periods for KB (don't have multiple observations to average)


#Explore results with plots --> IMPORT IN WITH NEW DATES
beta<- read.csv("euk_temporalBD_15s.csv") #bsaed on computations with RAW data
beta.long<- melt(beta, id.vars=c("Lake", "Comparison_type", "Comparison", "Comparison_years", "Comparison_years_CIC", "Year_midpoint_comparison", "Year_midpoint_comparisonCIC", "Comparison_order"))
colnames(beta.long)[9]<- 'Beta_component'
colnames(beta.long)[10]<- 'Prop'
as.factor(beta.long$Comparison_order)


#DAR
beta.long.DAR<- as.data.frame(subset(beta.long, Lake == "DAR", drop=T))
beta.long.DAR<- as.data.frame(subset(beta.long.DAR, Comparison_type == "Interval", drop=T))
#don't need non-existent time period comparisons 

betaDAR.plot<- ggplot(beta.long.DAR, aes(x=Comparison_order, y=Prop, colour = Beta_component)) + geom_point(size=4) + geom_path(size=2)
betaDAR.plot<- betaDAR.plot + theme_bw() + labs(x= "Comparison", y= "Beta")
betaDAR.plot<- betaDAR.plot + theme(axis.text.x = element_text(colour="black", size=16, angle=45))
betaDAR.plot<- betaDAR.plot + theme(axis.text.y = element_text(colour="black", size=16))
betaDAR.plot<- betaDAR.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
betaDAR.plot<- betaDAR.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
betaDAR.plot<- betaDAR.plot + scale_x_discrete(labels=c("1930-1945","1945-1950","1950-1963","1963-1973", "1973-1981", "1981-1994", "1994-2004", "2004-2011"))
betaDAR.plot<- betaDAR.plot + ggtitle("Dauriat")

#KB
beta.long.KB<- as.data.frame(subset(beta.long, Lake == "KB", drop=T))
beta.long.KB<- as.data.frame(subset(beta.long.KB, Comparison_type == "Interval", drop=T))

betaKB.plot<- ggplot(beta.long.KB, aes(x=Comparison_order, y=Prop, colour = Beta_component)) + geom_point(size=4) + geom_path(size=2)
betaKB.plot<- betaKB.plot + theme_bw() + labs(x= "Comparison", y= "Beta")
betaKB.plot<- betaKB.plot + theme(axis.text.x = element_text(colour="black", size=16, angle=45))
betaKB.plot<- betaKB.plot + theme(axis.text.y = element_text(colour="black", size=16))
betaKB.plot<- betaKB.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
betaKB.plot<- betaKB.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
betaKB.plot<- betaKB.plot + scale_x_discrete(labels=c("Pre1850-Pre1850","Pre1850-1866","1866-1928","1928-1982", "1982-2012"))
betaKB.plot<- betaKB.plot + ggtitle("Knob")

##Combo plots with midpoint years for comparisons
beta.long.red<- as.data.frame(subset(beta.long, Beta_component == "Total_beta" | Beta_component == "Species_loss"))
beta.long.red<- as.data.frame(subset(beta.long.red, Comparison_type == "Interval", drop=T))

#Manuscript figure 6
bd.plot2<- ggplot(beta.long.red, aes(x=Year_midpoint_comparisonCIC, y=Prop, shape=Beta_component)) + geom_point(size=4) + geom_path(size=1) + coord_flip() + facet_wrap(~Lake)
bd.plot2<- bd.plot2 + labs(x= "Midpoint of comparison", y= "Beta value")
bd.plot2<- bd.plot2 + theme_bw()
bd.plot2<- bd.plot2 + scale_shape_manual(values = c("Total_beta" = 17, "Species_loss" = 19))
bd.plot2<- bd.plot2 + theme(axis.text.x = element_text(colour="black", size=16))
bd.plot2<- bd.plot2 + theme(axis.text.y = element_text(colour="black", size=16))
bd.plot2<- bd.plot2 + theme(axis.title.x = element_text(size = rel(2), angle=00))
bd.plot2<- bd.plot2 + theme(axis.title.y = element_text(size = rel(2), angle=90))
bd.plot2<- bd.plot2 + annotate("rect", xmin=1939, xmax=1977, ymin=0, ymax=1, alpha=.2)
bd.plot2<- bd.plot2 + theme(strip.text.x = element_text(size=25, face="bold"), panel.background=element_rect(fill="white"))

#Manuscript figure 7 - adding in panel of cladoceran beta diversity from chapter 3. 
#(A) NGS beta diversity 
bd.plotA<- ggplot(beta.long.red, aes(x=Year_midpoint_comparisonCIC, y=Prop, shape=Beta_component)) + geom_point(size=4) + geom_path(size=1) + coord_flip() + facet_wrap(~Lake)
bd.plotA<- bd.plotA + labs(x= "Midpoint of comparison", y= "Beta value")
bd.plotA<- bd.plotA + theme_bw()
bd.plotA<- bd.plotA + scale_shape_manual(values = c("Total_beta" = 17, "Species_loss" = 19))
bd.plotA<- bd.plotA + theme(axis.text.x = element_text(colour="black", size=16))
bd.plotA<- bd.plotA + theme(axis.text.y = element_text(colour="black", size=16))
bd.plotA<- bd.plotA + theme(axis.title.x = element_text(size = rel(2), angle=00))
bd.plotA<- bd.plotA + theme(axis.title.y = element_text(size = rel(2), angle=90))
bd.plotA<- bd.plotA + annotate("rect", xmin=1939, xmax=1977, ymin=0, ymax=1, alpha=.2)
bd.plotA<- bd.plotA + theme(strip.text.x = element_text(size=25, face="bold"), panel.background=element_rect(fill="white"))
bd.plotA<- bd.plotA + theme(legend.position="none")

#(B) Cladoceran beta diversity - Ch. 3 
bdtemp.interval.res<- read.csv("schefferville_bdtemp_byinterval_res_clad.csv") #Cladoceran data 

interval.long<- melt(bdtemp.interval.res, id.vars=c("Lake", "Interval_comparison", "Comparison_order", "Year_comparisons", "Year_midpoint_comparison", "Year_comparisonsCIC", "Year_midpoint_comparisonCIC"))
colnames(interval.long)[8]<- 'Beta_component'
colnames(interval.long)[9]<- 'Value'

#Subset
interval.long.DAR<- as.data.frame(subset(interval.long, Lake == "DAR", drop=T))
interval.long.KB<- as.data.frame(subset(interval.long, Lake == "KB", drop=T))

interval.long.red<- as.data.frame(subset(interval.long, Beta_component == "Total_beta" | Beta_component == "Species_loss"))


bd.plotB<- ggplot(interval.long.red, aes(x=Year_midpoint_comparisonCIC, y=Value, shape=Beta_component)) + geom_point(size=4, colour="darkblue") + geom_path(size=1, colour="darkblue") + coord_flip() + facet_wrap(~Lake)
bd.plotB<- bd.plotB + labs(x= "Midpoint of comparison", y= "Beta value")
bd.plotB<- bd.plotB + theme_bw()
bd.plotB<- bd.plotB + scale_shape_manual(values = c("Total_beta" = 17, "Species_loss" = 19))
bd.plotB<- bd.plotB + theme(axis.text.x = element_text(colour="black", size=16))
bd.plotB<- bd.plotB + theme(axis.text.y = element_text(colour="black", size=16))
bd.plotB<- bd.plotB + theme(axis.title.x = element_text(size = rel(2), angle=00))
bd.plotB<- bd.plotB + theme(axis.title.y = element_text(size = rel(2), angle=90))
bd.plotB<- bd.plotB + annotate("rect", xmin=1939, xmax=1977, ymin=0, ymax=0.8, alpha=.2)
bd.plotB<- bd.plotB + theme(strip.text.x = element_text(size=25, face="bold"), panel.background=element_rect(fill="white"))
bd.plotB<- bd.plotB + theme(legend.position="none")
 

plotAB<- grid.arrange(bd.plotA, bd.plotB, nrow=2)


##Create a figure where data for both lakes on same scale- can add diatom data to this as well. 

bd.dat<- read.csv("Fig6_temporalBD_clad_andeuk.csv")

bd.dat.DAR<- as.data.frame(subset(bd.dat, Lake == "DAR", drop=T))
bd.dat.KB<- as.data.frame(subset(bd.dat, Lake == "KB", drop=T))


#Plot (A) - DAR
bd.dat.DAR.long<- melt(bd.dat.DAR, id.vars=c("Lake", "Data_type", "Year_midpoint_comparisonCIC"))
colnames(bd.dat.DAR.long)[4]<- 'Beta_component'
colnames(bd.dat.DAR.long)[5]<- 'Prop'
bd.dat.DAR.long.red<- as.data.frame(subset(bd.dat.DAR.long, Beta_component == "Total_beta" | Beta_component == "Species_loss"))

bd.plotA<- ggplot(bd.dat.DAR.long.red, aes(x=Year_midpoint_comparisonCIC, y=Prop, shape=Beta_component)) + geom_point(size=4, colour="darkblue") + geom_path(size=1, colour="darkblue") + coord_flip() + facet_wrap(~Data_type)
bd.plotA<- bd.plotA + labs(x= "Midpoint of comparison", y= "Beta value")
bd.plotA<- bd.plotA + theme_bw()
bd.plotA<- bd.plotA + scale_shape_manual(values = c("Total_beta" = 17, "Species_loss" = 19))
bd.plotA<- bd.plotA + theme(axis.text.x = element_text(colour="black", size=16))
bd.plotA<- bd.plotA + theme(axis.text.y = element_text(colour="black", size=16))
bd.plotA<- bd.plotA + theme(axis.title.x = element_text(size = rel(2), angle=00))
bd.plotA<- bd.plotA + theme(axis.title.y = element_text(size = rel(2), angle=90))
bd.plotA<- bd.plotA + annotate("rect", xmin=1939, xmax=1977, ymin=0, ymax=1, alpha=.2)
bd.plotA<- bd.plotA + theme(strip.text.x = element_text(size=21, face="bold"), panel.background=element_rect(fill="white"))
bd.plotA<- bd.plotA + theme(legend.position="none")

#Plot (B) - KB
bd.dat.KB.long<- melt(bd.dat.KB, id.vars=c("Lake", "Data_type", "Year_midpoint_comparisonCIC"))
colnames(bd.dat.KB.long)[4]<- 'Beta_component'
colnames(bd.dat.KB.long)[5]<- 'Prop'
bd.dat.KB.long.red<- as.data.frame(subset(bd.dat.KB.long, Beta_component == "Total_beta" | Beta_component == "Species_loss"))

bd.plotB<- ggplot(bd.dat.KB.long.red, aes(x=Year_midpoint_comparisonCIC, y=Prop, shape=Beta_component)) + geom_point(size=4, colour="darkred") + geom_path(size=1, colour="darkred") + coord_flip() + facet_wrap(~Data_type)
bd.plotB<- bd.plotB + labs(x= "Midpoint of comparison", y= "Beta value")
bd.plotB<- bd.plotB + theme_bw()
bd.plotB<- bd.plotB + scale_shape_manual(values = c("Total_beta" = 17, "Species_loss" = 19))
bd.plotB<- bd.plotB + theme(axis.text.x = element_text(colour="black", size=16))
bd.plotB<- bd.plotB + theme(axis.text.y = element_text(colour="black", size=16))
bd.plotB<- bd.plotB + theme(axis.title.x = element_text(size = rel(2), angle=00))
bd.plotB<- bd.plotB + theme(axis.title.y = element_text(size = rel(2), angle=90))
bd.plotB<- bd.plotB + annotate("rect", xmin=1939, xmax=1977, ymin=0, ymax=1, alpha=.2)
bd.plotB<- bd.plotB + theme(strip.text.x = element_text(size=21, face="bold"), panel.background=element_rect(fill="white"))
bd.plotB<- bd.plotB + theme(legend.position="none")

combo.plot<- grid.arrange(bd.plotA, bd.plotB, nrow=2)



#####ALSO INTERESTED IN USING EMBL3 GROUPINGS, BUT FROM RAW DATA- SO CAN LOOK AT
#GROUP SPECIFIC BETA DIVERSITY. 




#################PROCRUSTES TO COMPARE DYNAMICS BETWEEN CLUSTERED GROUPS##################
##Use matrix a (sample by OTUs) and the groups from hclust to complete ordinations
#for each hclust group, extract species scores and compare dynamics in the two 
#different groups. 

a.grp1<- as.data.frame(subset(a, clustgrp == "1", drop=T))
#Grp 1 includes:
#DAR1-2012
#DAR2- 1998
#DAR3- 1986
#DAR4- 1972
#DAR7- 1911 (outlier?)
#KB1- 2012
#KB2- 1982

a.grp2<- as.data.frame(subset(a, clustgrp == "2", drop=T))
#Grp 2 includes: 
#DAR5- 1963
#DAR6- 1947
#DAR8- 1866
#DAR9- Pre-1850
#KB3- 1928
#KB4- 1866
#KB5- Pre-1850
#KB6- Pre-1850

#PCAs - with Hellinger transformations 
#Group 1
grp1.pca<- rda(decostand(a.grp1[,5:887], method="hellinger"))
summary(ggrp1.pca)
plot(grp1.pca)

#SPECIES (OTU) scores 
grp1.sc<- scores(grp1.pca, choices=1:4, display="species", scaling=0)

#Group 2
grp2.pca<- rda(decostand(a.grp2[,5:887], method="hellinger"))
summary(ggrp2.pca)
plot(grp2.pca)

#SPECIES (OTU) scores 
grp2.sc<- scores(grp2.pca, choices=1:4, display="species", scaling=0)

#Procrustes overlay
protest(grp1.sc,grp2.sc) 
plot(procrustes(grp1.sc,grp2.sc))
proc.res<- as.data.frame(residuals(procrustes(grp1.sc, grp2.sc)))
colnames(proc.res)[1]<- 'Residual'
#write.csv(proc.res, "grp_procrustes_residuals.csv")



#################PIGMENTATION AND FUNCTIONAL GROUPS##################
##Match up OTUs and EMBL3 designations for OTUs which whether they are 
#pigmented (pi) or non-pigmented (npi) or unknown. 
#Also split into rough functional groups: pigmented/parasites/heterotrophs 

func.grp<- read.csv("euk_15s741r_OTUs_funcgrps.csv")
#reads in 3 extra blank columns, remove those. 
func.grp<- func.grp[,-4]
func.grp<- func.grp[,-4]
func.grp<- func.grp[,-4]

#Look at distribution across the groups (not considering sample or lake differences yet)
pi.plot<- qplot(Pigmentation, data=func.grp, geom="histogram")
pi.plot<- pi.plot + labs(x="Pigmentation", y="No. of OTUs") + theme_bw()

grp.plot<- qplot(Func_grp, data=func.grp, geom="histogram")
grp.plot<- grp.plot + labs(x="Functional group", y="No. of OTUs") + theme_bw()

#Bind the pigmentation and functional groups to matrix with lakes and 
#samples so can look at make-up per lake and per sample. 
a<-read.table("euk_nosing_15s741r.txt",h=T,row.names=1)

func.grp<- as.data.frame(cbind(func.grp, a[,1:15]))

#Make func.grp long and then combine with descriptive columns. 
func.grp.long<- melt(func.grp)
colnames(func.grp.long)[4]<- 'Sample_ID'
colnames(func.grp.long)[5]<- 'No_reads'
#Can work with long-format? 

pi.plot2<- ggplot(func.grp.long, aes(x=Sample_ID, y=No_reads, fill = Pigmentation)) + geom_histogram(stat = 'identity')
#Need to figure out this plot so all groups summed together. 

grp.plot2<- ggplot(func.grp.long, aes(x=Sample_ID, y=No_reads, fill = Func_grp)) + geom_histogram(stat = 'identity')



#################COMMON OTUs BETWEEN THE TWO LAKES##################
#Matrices
a<- read.table("euk_nosing_15s741r.txt",h=T) #OTU by sample matrix with rarefied reads
#columns 2:16 are the sediment samples. 

div <- read.csv("euk_nosing_15s741r_metrics_expanded.csv") #Exported metrics with lake info and approx dates added. 
#diversity metrics 

#Bind div to matrices a, b, and c and transform matrices. 
#OTUs - no. reads
a<- a[,1:16] #only need OTUs and samples
a<- melt(a) #make long-format
a<- dcast(a, variable ~ OTUs) #cast back into transposed matrix. 
a<- as.data.frame(cbind(div[,2:5], a[,2:884])) #bind with extra info
rownames(a)<- as.character(a[,2])

#Total # of OTUs between samples from both lakes: 883 OTUs

#Number of OTUs that have reads in DAR samples
a.DAR<- as.data.frame(subset(a, Lake == "DAR", drop=T))
summary(colSums(a.DAR[,5:887])>0) #585

#Number of OTUs that have reads in KB samples
a.KB<- as.data.frame(subset(a, Lake == "KB", drop=T))
summary(colSums(a.KB[,5:887])>0) #466

#Find number of OTUs in common 
#Keep only OTUs in each lake that have colSums greater than 0. 

a.DAR.unique<- a.DAR[,colSums(a.DAR[,5:887]) > 0] #585 OTUs, 588 columns  
a.KB.unique<- a.KB[,colSums(a.KB[,5:887]) > 0] #466 OTUs, 467 columns
#NB: not sure why less descriptive columns are retained in DAR vs. KB.  

#Check to see which are in common. 
#Format DAR
DAR.OTUlong<- as.data.frame(melt(a.DAR.unique, id.vars=c("Lake", "Approx_year", "Approx_year2")))
DAR.OTUlist<- as.data.frame(DAR.OTUlong$variable)
DAR.OTUlist<- as.data.frame(DAR.OTUlist[!duplicated(DAR.OTUlist[,1]),])
colnames(DAR.OTUlist)[1]<- 'DAR.OTUs'

#Format KB
KB.OTUlong<- as.data.frame(melt(a.KB.unique, id.vars=c("Sample_ID")))
KB.OTUlist<- as.data.frame(KB.OTUlong$variable)
KB.OTUlist<- as.data.frame(KB.OTUlist[!duplicated(KB.OTUlist[,1]),])
colnames(KB.OTUlist)[1]<- 'KB.OTUs'

#Match these lists 
match.check<- (DAR.OTUlist[,1] %in% KB.OTUlist[,1]) 
summary(match.check) #only 168 OTUs in common??


#################RUNNING ANALYSES BY FUNCTIONAL GROUPS##################
#OTUs by sample
a<- read.table("euk_nosing_15s741r.txt",h=T) #OTU by sample matrix with rarefied reads

#functional groups
func.grp<- read.csv("euk_15s741r_OTUs_funcgrps.csv")
#reads in 3 extra blank columns, remove those. 
func.grp<- func.grp[,-4]
func.grp<- func.grp[,-4]
func.grp<- func.grp[,-4]

#diversity metrics with descriptive columns 
div <- read.csv("euk_nosing_15s741r_metrics_expanded.csv") #Exported metrics with lake info and approx dates added. 


#Subset func.grp into the 4 different functional groups 
func.grp.par<- as.data.frame(subset(func.grp, Func_grp == "Parasite", drop=T))
func.grp.het<- as.data.frame(subset(func.grp, Func_grp == "Het_other", drop=T))
func.grp.unk<- as.data.frame(subset(func.grp, Func_grp == "UNK", drop=T))
func.grp.auto<- as.data.frame(subset(func.grp, Func_grp == "Autotroph", drop=T))


#Match 'a' matrix each subseted func.grp.XX dataframe. 

#a.par
#Sites in 'a' that are found in func.grp.par
par.check<- (a$OTUs %in% func.grp.par$OTUs) 
a.par<- a[a$OTUs %in% func.grp.par$OTUs, ]
  
#a.het
#Sites in 'a' that are found in func.grp.het
a.het<- a[a$OTUs %in% func.grp.het$OTUs, ]

#a.unk
#Sites in 'a' that are found in func.grp.unk
a.unk<- a[a$OTUs %in% func.grp.unk$OTUs, ]

#a.auto
#Sites in 'a' that are found in func.grp.auto
a.auto<- a[a$OTUs %in% func.grp.auto$OTUs, ]


#Also create matrices of these OTUs that are sample by OTU with additional columns for 
#descriptive variables. 

#a.par
a.par<- a.par[,1:16] 
a.par<- melt(a.par) 
a.par<- dcast(a.par, variable ~ OTUs)
a.par<- as.data.frame(cbind(div[,2:5], a.par[,2:140])) 
rownames(a.par)<- as.character(a.par[,2])

#a.het
a.het<- a.het[,1:16] 
a.het<- melt(a.het) 
a.het<- dcast(a.het, variable ~ OTUs)
a.het<- as.data.frame(cbind(div[,2:5], a.het[,2:75])) 
rownames(a.het)<- as.character(a.het[,2])

#a.unk
a.unk<- a.unk[,1:16] 
a.unk<- melt(a.unk) 
a.unk<- dcast(a.unk, variable ~ OTUs)
a.unk<- as.data.frame(cbind(div[,2:5], a.unk[,2:370])) 
rownames(a.unk)<- as.character(a.unk[,2])

#a.auto 
a.auto<- a.auto[,1:16] 
a.auto<- melt(a.auto) 
a.auto<- dcast(a.auto, variable ~ OTUs)
a.auto<- as.data.frame(cbind(div[,2:5], a.auto[,2:302])) 
rownames(a.auto)<- as.character(a.auto[,2])

#Matrices to use for analyses:
a.par #parasites, 139 OTUs
a.het #heterotrophs, 74 OTUs
a.unk #unknown functional group, 369 OTUs
a.auto #autotrophs , 301 OTUs 

##ANALYSES WITH THESE FOUR DIFFERENT GROUPS 

####
##PARASITES
#NMDS
#NMDS - sample by OTU (No. reads)
par.nmds<-metaMDS(a.par[,5:143],distance="bray") 
par1 <- par.nmds$points[,1]
par2 <- par.nmds$points[,2]
parc<- as.data.frame(cbind(par1, par2))

par.nmds.plot<- ggplot()
par.nmds.plot<- par.nmds.plot + geom_vline(x=0,colour="grey50")
par.nmds.plot<- par.nmds.plot + geom_hline(y=0,colour="grey50")
par.nmds.plot<- par.nmds.plot + geom_text(data = parc, aes(x=par1, y=par2, label = a.par$Approx_year2, colour = a.par$Lake), size = 7)
par.nmds.plot<- par.nmds.plot + scale_colour_manual(values = c("DAR" = "darkblue", "KB" = "red2"), "Lake")
par.nmds.plot<- par.nmds.plot + theme_bw()
par.nmds.plot<- par.nmds.plot + labs(x= "NMDS 1", y= "NMDS 2")
par.nmds.plot<- par.nmds.plot + theme(axis.text.x = element_text(colour="black", size=16))
par.nmds.plot<- par.nmds.plot + theme(axis.text.y = element_text(colour="black", size=16))
par.nmds.plot<- par.nmds.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
par.nmds.plot<- par.nmds.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
par.nmds.plot<- par.nmds.plot + ggtitle("Parasite NMDS - OTUs (Reads non-transformed)")
#Add if want vectors
par.nmds.plot<- par.nmds.plot + geom_path(data = parc, aes(x = par1, y= par2, colour= a.par$Lake), size =0.4, arrow=arrow(length= unit(0.4,"cm"), ends="first", type="open" ))

#NMDS with hclust
par.bray<-vegdist(a.par[,5:143], method="bray")
par.clust<-hclust(par.bray)
par.grp<- cutree(par.clust, 3) #cutting into 3 groups - can see 3 in dendogram

par.nmds<-metaMDS(a.par[,5:143],distance="bray") 

col<- c("red2", "mediumblue", "black")
col[par.grp]

plot(par.nmds, type = "n", display = "sites")
points(par.nmds, col = col[par.grp], bg = col[par.grp], pch = 21)
legend("topright", legend = paste("Cluster", 1:2),
       col = col, pt.bg = col, bty = "n", pch = 21)

plot(par.nmds, type = "n", display = "sites")
text(par.nmds, col = col[par.grp], bg = col[par.grp], label=a.par$Sample_ID)
legend("topright", legend = paste("Cluster", 1:2),
       col = col, pt.bg = col, bty = "n", pch = 21)

par1 <- par.nmds$points[,1]
par2 <- par.nmds$points[,2]
parc<- as.data.frame(cbind(par1, par2, par.grp))
colnames(parc)[3]<- 'grp'
pgrp<- as.factor(parc$grp)

par.nmds.plot<- ggplot()
par.nmds.plot<- par.nmds.plot + geom_vline(x=0,colour="grey50")
par.nmds.plot<- par.nmds.plot + geom_hline(y=0,colour="grey50")
par.nmds.plot<- par.nmds.plot + geom_text(data = parc, aes(x=par1, y=par2, label = a.par$Approx_year2, colour = pgrp), size = 7)
par.nmds.plot<- par.nmds.plot + theme_bw()
par.nmds.plot<- par.nmds.plot + labs(x= "NMDS 1", y= "NMDS 2")
par.nmds.plot<- par.nmds.plot + theme(axis.text.x = element_text(colour="black", size=16))
par.nmds.plot<- par.nmds.plot + theme(axis.text.y = element_text(colour="black", size=16))
par.nmds.plot<- par.nmds.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
par.nmds.plot<- par.nmds.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
par.nmds.plot<- par.nmds.plot + ggtitle("Parasite NMDS - OTUs (Reads non-transformed, with clusters)")
#Add if want vectors
par.nmds.plot<- par.nmds.plot + geom_path(data = parc, aes(x = par1, y= par2, line = a.par$Lake), size =0.4, arrow=arrow(length= unit(0.4,"cm"), ends="first", type="open" ))

#ANOSIM
#Between lakes
par.ano<- anosim(a.par[,5:143], a.par$Lake, permutations = 999, distance = "bray")
summary(par.ano)

#Between Time periods 
a.par$Time_period<- c("Post-mining", "Post-mining", "Post-mining", "Mining", "Mining", "Mining", "Pre-mining", "Pre-mining", "Pre-mining", "Post-mining", "Mining", "Post-mining", "Post-mining", "Post-mining", "Post-mining")
par.ano<- anosim(a.par[,5:143], a.par$Time_period, permutations = 999, distance = "bray")
summary(par.ano)

#Between hclust groups - added manually from hclust results 
a.par$hclust<- c(1,1,1,1,1,2,1,1,1,3,3,2,2,3,2)
par.ano<- anosim(a.par[,5:143], a.par$hclust, permutations = 999, distance = "bray")
summary(par.ano)

#MRTs
#Metal data
met<- read.table("euk_metalEF.csv",h=T, sep=",")
a.par$met<- met$Metal_EF

#Year MRT
#Subset into the two lakes
a.par.DAR<- as.data.frame(subset(a.par, Lake == 'DAR', drop=T))
a.par.KB<- as.data.frame(subset(a.par, Lake == 'KB', drop=T))

y<- as.matrix(a.par.KB[,5:143]) #swap out a.DAR with a.KB here and in mvpart() (next line)

par.mrt<- mvpart(y ~ a.par.KB$Approx_year2, data = a.par.KB, cp=0, xv="pick", margin=0.08, xvmult=100, which=4) 
par.mrt
printcp(par.mrt) #Approx_year2 used
plotcp(par.mrt,minline=TRUE,upper=c("size"))
summary(par.mrt)
rsq.rpart(par.mrt)
par.mrt$cptable
1-par.mrt$cptable[2,3] 
#Plot
tree.par.mrt <- MRT(par.mrt,percent=10,species=colnames(y))
plot(tree.par.mrt)
summary(tree.par.mrt)

#Metal MRT
y<- as.matrix(a.par[,5:143]) 

par.mrt<- mvpart(y ~ a.par$met, data = a.par, cp=0, xv="pick", margin=0.08, xvmult=100, which=4) 
par.mrt
printcp(par.mrt) #metal used
plotcp(par.mrt,minline=TRUE,upper=c("size"))
summary(par.mrt)
rsq.rpart(par.mrt)
par.mrt$cptable
1-par.mrt$cptable[4,3] 
#Plot
tree.par.mrt <- MRT(par.mrt,percent=10,species=colnames(y))
plot(tree.par.mrt)
summary(tree.par.mrt)

#SIMPER - do in PAST


####
##HETEROTROPHS
#NMDS
#NMDS - sample by OTU (No. reads)
het.nmds<-metaMDS(a.het[,5:78],distance="bray") 
het1 <- het.nmds$points[,1]
het2 <- het.nmds$points[,2]
hetc<- as.data.frame(cbind(het1, het2))

het.nmds.plot<- ggplot()
het.nmds.plot<- het.nmds.plot + geom_vline(x=0,colour="grey50")
het.nmds.plot<- het.nmds.plot + geom_hline(y=0,colour="grey50")
het.nmds.plot<- het.nmds.plot + geom_text(data = hetc, aes(x=het1, y=het2, label = a.het$Approx_year2, colour = a.het$Lake), size = 7)
het.nmds.plot<- het.nmds.plot + scale_colour_manual(values = c("DAR" = "darkblue", "KB" = "red2"), "Lake")
het.nmds.plot<- het.nmds.plot + theme_bw()
het.nmds.plot<- het.nmds.plot + labs(x= "NMDS 1", y= "NMDS 2")
het.nmds.plot<- het.nmds.plot + theme(axis.text.x = element_text(colour="black", size=16))
het.nmds.plot<- het.nmds.plot + theme(axis.text.y = element_text(colour="black", size=16))
het.nmds.plot<- het.nmds.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
het.nmds.plot<- het.nmds.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
het.nmds.plot<- het.nmds.plot + ggtitle("Heterotroph NMDS - OTUs (Reads non-transformed)")
#Add if want vectors
het.nmds.plot<- het.nmds.plot + geom_path(data = hetc, aes(x = het1, y= het2, colour= a.het$Lake), size =0.4, arrow=arrow(length= unit(0.4,"cm"), ends="first", type="open" ))

#NMDS with hclust
het.bray<-vegdist(a.par[,5:78], method="bray")
het.clust<-hclust(het.bray)
het.grp<- cutree(het.clust, 3) #cutting into 3 groups - can see 3 in dendogram

het.nmds<-metaMDS(a.het[,5:78],distance="bray") 

col<- c("red2", "mediumblue", "black")
col[het.grp]

plot(het.nmds, type = "n", display = "sites")
points(het.nmds, col = col[het.grp], bg = col[het.grp], pch = 21)
legend("topright", legend = paste("Cluster", 1:2),
       col = col, pt.bg = col, bty = "n", pch = 21)

plot(het.nmds, type = "n", display = "sites")
text(het.nmds, col = col[het.grp], bg = col[het.grp], label=a.het$Sample_ID)
legend("topright", legend = paste("Cluster", 1:2),
       col = col, pt.bg = col, bty = "n", pch = 21)

het1 <- het.nmds$points[,1]
het2 <- het.nmds$points[,2]
hetc<- as.data.frame(cbind(het1, het2, het.grp))
colnames(hetc)[3]<- 'grp'
hgrp<- as.factor(hetc$grp)

het.nmds.plot<- ggplot()
het.nmds.plot<- het.nmds.plot + geom_vline(x=0,colour="grey50")
het.nmds.plot<- het.nmds.plot + geom_hline(y=0,colour="grey50")
het.nmds.plot<- het.nmds.plot + geom_text(data = hetc, aes(x=het1, y=het2, label = a.het$Approx_year2, colour = hgrp), size = 7)
het.nmds.plot<- het.nmds.plot + theme_bw()
het.nmds.plot<- het.nmds.plot + labs(x= "NMDS 1", y= "NMDS 2")
het.nmds.plot<- het.nmds.plot + theme(axis.text.x = element_text(colour="black", size=16))
het.nmds.plot<- het.nmds.plot + theme(axis.text.y = element_text(colour="black", size=16))
het.nmds.plot<- het.nmds.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
het.nmds.plot<- het.nmds.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
het.nmds.plot<- het.nmds.plot + ggtitle("Heterotroph NMDS - OTUs (Reads non-transformed, with clusters)")
#Add if want vectors
het.nmds.plot<- het.nmds.plot + geom_path(data = hetc, aes(x = het1, y= het2, line = a.het$Lake), size =0.4, arrow=arrow(length= unit(0.4,"cm"), ends="first", type="open" ))

#ANOSIM
#Between lakes
het.ano<- anosim(a.het[,5:78], a.het$Lake, permutations = 999, distance = "bray")
summary(het.ano)

#Between Time periods 
a.het$Time_period<- c("Post-mining", "Post-mining", "Post-mining", "Mining", "Mining", "Mining", "Pre-mining", "Pre-mining", "Pre-mining", "Post-mining", "Mining", "Post-mining", "Post-mining", "Post-mining", "Post-mining")
het.ano<- anosim(a.het[,5:78], a.het$Time_period, permutations = 999, distance = "bray")
summary(het.ano)

#Between hclust groups - added manually from hclust results 
a.het$hclust<- c(1,1,2,2,2,2,1,1,2,3,3,2,2,3,2)
het.ano<- anosim(a.het[,5:78], a.het$hclust, permutations = 999, distance = "bray")
summary(het.ano)

#MRTs
#Metal data
a.het$met<- met$Metal_EF

#Year MRT
#Subset into the two lakes
a.het.DAR<- as.data.frame(subset(a.het, Lake == 'DAR', drop=T))
a.het.KB<- as.data.frame(subset(a.het, Lake == 'KB', drop=T))

y<- as.matrix(a.het.KB[,5:78]) #swap out a.DAR with a.KB here and in mvpart() (next line)

het.mrt<- mvpart(y ~ a.het.KB$Approx_year2, data = a.het.KB, cp=0, xv="pick", margin=0.08, xvmult=100, which=4) 
het.mrt
printcp(het.mrt) #Approx_year2 used
plotcp(het.mrt,minline=TRUE,upper=c("size"))
summary(het.mrt)
rsq.rpart(het.mrt)
het.mrt$cptable
1-het.mrt$cptable[2,3] 
#Plot
tree.het.mrt <- MRT(het.mrt,percent=10,species=colnames(y))
plot(tree.het.mrt)
summary(tree.het.mrt)

#Metal MRT
y<- as.matrix(a.het[,5:78]) 

het.mrt<- mvpart(y ~ a.het$met, data = a.het, cp=0, xv="pick", margin=0.08, xvmult=100, which=4) 
het.mrt
printcp(het.mrt) #metal used
plotcp(het.mrt,minline=TRUE,upper=c("size"))
summary(het.mrt)
rsq.rpart(het.mrt)
het.mrt$cptable
1-het.mrt$cptable[4,3] 
#Plot
tree.het.mrt <- MRT(het.mrt,percent=10,species=colnames(y))
plot(tree.het.mrt)
summary(tree.het.mrt)

#SIMPER - do in PAST


####
##UNKNOWN
#NMDS
#NMDS - sample by OTU (No. reads)
unk.nmds<-metaMDS(a.unk[,5:373],distance="bray") 
unk1 <- unk.nmds$points[,1]
unk2 <- unk.nmds$points[,2]
unkc<- as.data.frame(cbind(unk1, unk2))

unk.nmds.plot<- ggplot()
unk.nmds.plot<- unk.nmds.plot + geom_vline(x=0,colour="grey50")
unk.nmds.plot<- unk.nmds.plot + geom_hline(y=0,colour="grey50")
unk.nmds.plot<- unk.nmds.plot + geom_text(data = unkc, aes(x=unk1, y=unk2, label = a.unk$Approx_year2, colour = a.unk$Lake), size = 7)
unk.nmds.plot<- unk.nmds.plot + scale_colour_manual(values = c("DAR" = "darkblue", "KB" = "red2"), "Lake")
unk.nmds.plot<- unk.nmds.plot + theme_bw()
unk.nmds.plot<- unk.nmds.plot + labs(x= "NMDS 1", y= "NMDS 2")
unk.nmds.plot<- unk.nmds.plot + theme(axis.text.x = element_text(colour="black", size=16))
unk.nmds.plot<- unk.nmds.plot + theme(axis.text.y = element_text(colour="black", size=16))
unk.nmds.plot<- unk.nmds.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
unk.nmds.plot<- unk.nmds.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
unk.nmds.plot<- unk.nmds.plot + ggtitle("Unclassified NMDS - OTUs (Reads non-transformed)")
#Add if want vectors
unk.nmds.plot<- unk.nmds.plot + geom_path(data = unkc, aes(x = unk1, y= unk2, colour= a.unk$Lake), size =0.4, arrow=arrow(length= unit(0.4,"cm"), ends="first", type="open" ))

#NMDS with hclust
unk.bray<-vegdist(a.unk[,5:373], method="bray")
unk.clust<-hclust(unk.bray)
unk.grp<- cutree(unk.clust, 2) #cutting into 2 groups - can see 2 in dendogram

unk.nmds<-metaMDS(a.unk[,5:373],distance="bray") 

col<- c("red2", "mediumblue", "black")
col[unk.grp]

plot(unk.nmds, type = "n", display = "sites")
points(unk.nmds, col = col[unk.grp], bg = col[unk.grp], pch = 21)
legend("topright", legend = paste("Cluster", 1:2),
       col = col, pt.bg = col, bty = "n", pch = 21)

plot(unk.nmds, type = "n", display = "sites")
text(unk.nmds, col = col[unk.grp], bg = col[unk.grp], label=a.unk$Sample_ID)
legend("topright", legend = paste("Cluster", 1:2),
       col = col, pt.bg = col, bty = "n", pch = 21)

unk1 <- unk.nmds$points[,1]
unk2 <- unk.nmds$points[,2]
unkc<- as.data.frame(cbind(unk1, unk2, unk.grp))
colnames(unkc)[3]<- 'grp'
ugrp<- as.factor(unkc$grp)

unk.nmds.plot<- ggplot()
unk.nmds.plot<- unk.nmds.plot + geom_vline(x=0,colour="grey50")
unk.nmds.plot<- unk.nmds.plot + geom_hline(y=0,colour="grey50")
unk.nmds.plot<- unk.nmds.plot + geom_text(data = unkc, aes(x=unk1, y=unk2, label = a.unk$Approx_year2, colour = ugrp), size = 7)
unk.nmds.plot<- unk.nmds.plot + theme_bw()
unk.nmds.plot<- unk.nmds.plot + labs(x= "NMDS 1", y= "NMDS 2")
unk.nmds.plot<- unk.nmds.plot + theme(axis.text.x = element_text(colour="black", size=16))
unk.nmds.plot<- unk.nmds.plot + theme(axis.text.y = element_text(colour="black", size=16))
unk.nmds.plot<- unk.nmds.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
unk.nmds.plot<- unk.nmds.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
unk.nmds.plot<- unk.nmds.plot + ggtitle("Unclassified NMDS - OTUs (Reads non-transformed, with clusters)")
#Add if want vectors
unk.nmds.plot<- unk.nmds.plot + geom_path(data = unkc, aes(x = unk1, y= unk2, line = a.unk$Lake), size =0.4, arrow=arrow(length= unit(0.4,"cm"), ends="first", type="open" ))

#ANOSIM
#Between lakes
unk.ano<- anosim(a.unk[,5:373], a.unk$Lake, permutations = 999, distance = "bray")
summary(unk.ano)

#Between Time periods 
a.unk$Time_period<- c("Post-mining", "Post-mining", "Post-mining", "Mining", "Mining", "Mining", "Pre-mining", "Pre-mining", "Pre-mining", "Post-mining", "Mining", "Post-mining", "Post-mining", "Post-mining", "Post-mining")
unk.ano<- anosim(a.unk[,5:373], a.unk$Time_period, permutations = 999, distance = "bray")
summary(unk.ano)

#Between hclust groups - added manually from hclust results 
a.unk$hclust<- c(1,1,1,1,1,1,1,1,1,2,2,2,2,2,2)
unk.ano<- anosim(a.unk[,5:373], a.unk$hclust, permutations = 999, distance = "bray")
summary(unk.ano)

#MRTs
#Metal data
a.unk$met<- met$Metal_EF

#Year MRT
#Subset into the two lakes
a.unk.DAR<- as.data.frame(subset(a.unk, Lake == 'DAR', drop=T))
a.unk.KB<- as.data.frame(subset(a.unk, Lake == 'KB', drop=T))

y<- as.matrix(a.unk.KB[,5:373]) #swap out a.DAR with a.KB here and in mvpart() (next line)

unk.mrt<- mvpart(y ~ a.unk.KB$Approx_year2, data = a.unk.KB, cp=0, xv="pick", margin=0.08, xvmult=100, which=4) 
unk.mrt
printcp(unk.mrt) #Approx_year2 used
plotcp(unk.mrt,minline=TRUE,upper=c("size"))
summary(unk.mrt)
rsq.rpart(unk.mrt)
unk.mrt$cptable
1-unk.mrt$cptable[2,3] 
#Plot
tree.unk.mrt <- MRT(unk.mrt,percent=10,species=colnames(y))
plot(tree.unk.mrt)
summary(tree.unk.mrt)

#Metal MRT
y<- as.matrix(a.unk[,5:373]) 

unk.mrt<- mvpart(y ~ a.unk$met, data = a.unk, cp=0, xv="pick", margin=0.08, xvmult=100, which=4) 
unk.mrt
printcp(unk.mrt) #metal used
plotcp(unk.mrt,minline=TRUE,upper=c("size"))
summary(unk.mrt)
rsq.rpart(unk.mrt)
unk.mrt$cptable
1-unk.mrt$cptable[4,3] 
#Plot
tree.unk.mrt <- MRT(unk.mrt,percent=10,species=colnames(y))
plot(tree.unk.mrt)
summary(tree.unk.mrt)

#SIMPER - do in PAST


####
##AUTOTROPHS
#NMDS
#NMDS - sample by OTU (No. reads)
auto.nmds<-metaMDS(a.auto[,5:305],distance="bray") 
auto1 <- auto.nmds$points[,1]
auto2 <- auto.nmds$points[,2]
autoc<- as.data.frame(cbind(auto1, auto2))

auto.nmds.plot<- ggplot()
auto.nmds.plot<- auto.nmds.plot + geom_vline(x=0,colour="grey50")
auto.nmds.plot<- auto.nmds.plot + geom_hline(y=0,colour="grey50")
auto.nmds.plot<- auto.nmds.plot + geom_text(data = autoc, aes(x=auto1, y=auto2, label = a.auto$Approx_year2, colour = a.auto$Lake), size = 7)
auto.nmds.plot<- auto.nmds.plot + scale_colour_manual(values = c("DAR" = "darkblue", "KB" = "red2"), "Lake")
auto.nmds.plot<- auto.nmds.plot + theme_bw()
auto.nmds.plot<- auto.nmds.plot + labs(x= "NMDS 1", y= "NMDS 2")
auto.nmds.plot<- auto.nmds.plot + theme(axis.text.x = element_text(colour="black", size=16))
auto.nmds.plot<- auto.nmds.plot + theme(axis.text.y = element_text(colour="black", size=16))
auto.nmds.plot<- auto.nmds.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
auto.nmds.plot<- auto.nmds.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
auto.nmds.plot<- auto.nmds.plot + ggtitle("Autotroph NMDS - OTUs (Reads non-transformed)")
#Add if want vectors
auto.nmds.plot<- auto.nmds.plot + geom_path(data = autoc, aes(x = auto1, y= auto2, colour= a.auto$Lake), size =0.4, arrow=arrow(length= unit(0.4,"cm"), ends="first", type="open" ))

#NMDS with hclust
auto.bray<-vegdist(a.auto[,5:305], method="bray")
auto.clust<-hclust(auto.bray)
auto.grp<- cutree(auto.clust, 3) #cutting into 3 groups - can see 3 in dendogram

auto.nmds<-metaMDS(a.auto[,5:305],distance="bray") 

col<- c("red2", "mediumblue", "black")
col[auto.grp]

plot(auto.nmds, type = "n", display = "sites")
points(auto.nmds, col = col[auto.grp], bg = col[auto.grp], pch = 21)
legend("topright", legend = paste("Cluster", 1:2),
       col = col, pt.bg = col, bty = "n", pch = 21)

plot(auto.nmds, type = "n", display = "sites")
text(auto.nmds, col = col[auto.grp], bg = col[auto.grp], label=a.auto$Sample_ID)
legend("topright", legend = paste("Cluster", 1:2),
       col = col, pt.bg = col, bty = "n", pch = 21)

auto1 <- auto.nmds$points[,1]
auto2 <- auto.nmds$points[,2]
autoc<- as.data.frame(cbind(auto1, auto2, auto.grp))
colnames(autoc)[3]<- 'grp'
agrp<- as.factor(autoc$grp)

auto.nmds.plot<- ggplot()
auto.nmds.plot<- auto.nmds.plot + geom_vline(x=0,colour="grey50")
auto.nmds.plot<- auto.nmds.plot + geom_hline(y=0,colour="grey50")
auto.nmds.plot<- auto.nmds.plot + geom_text(data = autoc, aes(x=auto1, y=auto2, label = a.auto$Approx_year2, colour = agrp), size = 7)
auto.nmds.plot<- auto.nmds.plot + theme_bw()
auto.nmds.plot<- auto.nmds.plot + labs(x= "NMDS 1", y= "NMDS 2")
auto.nmds.plot<- auto.nmds.plot + theme(axis.text.x = element_text(colour="black", size=16))
auto.nmds.plot<- auto.nmds.plot + theme(axis.text.y = element_text(colour="black", size=16))
auto.nmds.plot<- auto.nmds.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
auto.nmds.plot<- auto.nmds.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
auto.nmds.plot<- auto.nmds.plot + ggtitle("Autotroph NMDS - OTUs (Reads non-transformed, with clusters)")
#Add if want vectors
auto.nmds.plot<- auto.nmds.plot + geom_path(data = autoc, aes(x = auto1, y= auto2, line = a.auto$Lake), size =0.4, arrow=arrow(length= unit(0.4,"cm"), ends="first", type="open" ))

#ANOSIM
#Between lakes
auto.ano<- anosim(a.auto[,5:305], a.auto$Lake, permutations = 999, distance = "bray")
summary(auto.ano)

#Between Time periods 
a.auto$Time_period<- c("Post-mining", "Post-mining", "Post-mining", "Mining", "Mining", "Mining", "Pre-mining", "Pre-mining", "Pre-mining", "Post-mining", "Mining", "Post-mining", "Post-mining", "Post-mining", "Post-mining")
auto.ano<- anosim(a.auto[,5:305], a.auto$Time_period, permutations = 999, distance = "bray")
summary(auto.ano)

#Between hclust groups - added manually from hclust results 
a.auto$hclust<- c(1,1,2,1,1,2,1,2,3,1,1,2,2,2,2)
auto.ano<- anosim(a.auto[,5:305], a.auto$hclust, permutations = 999, distance = "bray")
summary(auto.ano)

#MRTs
#Metal data
a.auto$met<- met$Metal_EF

#Year MRT
#Subset into the two lakes
a.auto.DAR<- as.data.frame(subset(a.auto, Lake == 'DAR', drop=T))
a.auto.KB<- as.data.frame(subset(a.auto, Lake == 'KB', drop=T))

y<- as.matrix(a.auto.KB[,5:305]) #swap out a.DAR with a.KB here and in mvpart() (next line)

auto.mrt<- mvpart(y ~ a.auto.KB$Approx_year2, data = a.auto.KB, cp=0, xv="pick", margin=0.08, xvmult=100, which=4) 
auto.mrt
printcp(auto.mrt) #Approx_year2 used
plotcp(auto.mrt,minline=TRUE,upper=c("size"))
summary(auto.mrt)
rsq.rpart(auto.mrt)
auto.mrt$cptable
1-auto.mrt$cptable[2,3] 
#Plot
tree.auto.mrt <- MRT(auto.mrt,percent=10,species=colnames(y))
plot(tree.auto.mrt)
summary(tree.auto.mrt)

#Metal MRT
y<- as.matrix(a.auto[,5:305]) 

auto.mrt<- mvpart(y ~ a.auto$met, data = a.auto, cp=0, xv="pick", margin=0.08, xvmult=100, which=4) 
auto.mrt
printcp(auto.mrt) #metal used
plotcp(auto.mrt,minline=TRUE,upper=c("size"))
summary(auto.mrt)
rsq.rpart(auto.mrt)
auto.mrt$cptable
1-auto.mrt$cptable[4,3] 
#Plot
tree.auto.mrt <- MRT(auto.mrt,percent=10,species=colnames(y))
plot(tree.auto.mrt)
summary(tree.auto.mrt)

#SIMPER - do in PAST


#################CHECKING RESULTS USING VALVE DIATOM DATA (LAPERRIERE)#################
##MRTs
##Using diatom valve data from Laperriere et al. 
#In order to check breakpoints uncovered by MRTs and to complete on dataset with continuous years. 

#Diatom valve data (counts)
diat<- read.csv("dauriat_subfossil_diatoms1998.csv") 
  
#Rarefy counts
raremax<-min(rowSums(diat[,6:62]))
Srare<-as.data.frame(rarefy(diat[,6:62],raremax)) #ok, may use later on. 

#Hellinger transform
diat.hell<- as.data.frame(decostand(diat[,6:62], method = "hell"))
diat.hell<- as.data.frame(cbind(diat[,1:5], diat.hell))

#MRTs
y <- as.matrix(diat.hell[,6:62])
a.mrt<- mvpart(y ~ diat.hell$Estimated_year, data = diat.hell, cp=0, xv="pick", margin=0.08, xvmult=100, which=4) 
a.mrt
a.mrt$frame$dev
printcp(a.mrt) 
plotcp(a.mrt,minline=TRUE,upper=c("size"))
summary(a.mrt)
par(mfrow=c(1,2)) # two plots on one page 
rsq.rpart(a.mrt)  # visualize cross-validation result
#R2 is 1-Rel Error
# R2  
a.mrt$cptable
1-a.mrt$cptable[12,3] #0.7

#Plot (MVPARTwrap pkg)
par(mfrow=c(1,1))
tree.a.mrt <- MRT(a.mrt,percent=10,species=colnames(y))
plot(tree.a.mrt)
summary(tree.a.mrt)


##NMDS and HCLUST using diatom valve data.
#Cut into 2 groups
a.bray<-vegdist(diat[,6:62], method="bray")
a.clust<-hclust(a.bray)
a.grp<- cutree(a.clust, 2) #cutting into 2 groups 

a.nmds<-metaMDS(diat[,6:62],distance="bray") 

col<- c("red2", "mediumblue")
col[a.grp]

plot(a.nmds, type = "n", display = "sites")
points(a.nmds, col = col[a.grp], bg = col[a.grp], pch = 21)
legend("topright", legend = paste("Cluster", 1:2),
       col = col, pt.bg = col, bty = "n", pch = 21)

plot(a.nmds, type = "n", display = "sites")
text(a.nmds, col = col[a.grp], bg = col[a.grp], label=diat$Estimated_year)
legend("topright", legend = paste("Cluster", 1:2),
       col = col, pt.bg = col, bty = "n", pch = 21)

#Try with ggplot so can better visualize 
a.bray<-vegdist(diat[,6:62], method="bray") #distance matrix 
a.clust<-hclust(a.bray) #cluster
a.grp<- as.data.frame(cutree(a.clust, 2)) #define groups 
a.nmds<-metaMDS(diat[,6:62],distance="bray") 
a1 <- a.nmds$points[,1]
a2 <- a.nmds$points[,2]
ac<- as.data.frame(cbind(a1, a2, a.grp))
colnames(ac)[3]<- 'grp'
grp<- as.factor(ac$grp)

a.nmds.plot<- ggplot()
a.nmds.plot<- a.nmds.plot + geom_vline(x=0,colour="grey50")
a.nmds.plot<- a.nmds.plot + geom_hline(y=0,colour="grey50")
a.nmds.plot<- a.nmds.plot + geom_text(data = ac, aes(x=a1, y=a2, label = diat$Estimated_year, colour = grp), size = 7)
#a.nmds.plot<- a.nmds.plot + scale_colour_manual(values = c("1" = "darkblue", "2" = "red2"), "Cluster group")
a.nmds.plot<- a.nmds.plot + theme_bw()
a.nmds.plot<- a.nmds.plot + labs(x= "NMDS 1", y= "NMDS 2")
a.nmds.plot<- a.nmds.plot + theme(axis.text.x = element_text(colour="black", size=16))
a.nmds.plot<- a.nmds.plot + theme(axis.text.y = element_text(colour="black", size=16))
a.nmds.plot<- a.nmds.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
a.nmds.plot<- a.nmds.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
#a.nmds.plot<- a.nmds.plot + ggtitle("NMDS - Counts")


#Cutting into 3 groups 
b.bray<-vegdist(diat[,6:62], method="bray")
b.clust<-hclust(b.bray)
b.grp<- cutree(b.clust, 3) #cutting into 3 groups 

b.nmds<-metaMDS(diat[,6:62],distance="bray") 

col<- c("red2", "mediumblue", "black")
col[b.grp]

plot(b.nmds, type = "n", display = "sites")
points(b.nmds, col = col[b.grp], bg = col[b.grp], pch = 21)
legend("topright", legend = paste("Cluster", 1:2),
       col = col, pt.bg = col, bty = "n", pch = 21)

plot(b.nmds, type = "n", display = "sites")
text(b.nmds, col = col[b.grp], bg = col[b.grp], label=diat$Estimated_year)
legend("topright", legend = paste("Cluster", 1:3),
       col = col, pt.bg = col, bty = "n", pch = 21)

#Try with ggplot so can better visualize 
b.bray<-vegdist(diat[,6:62], method="bray") #distance matrix 
b.clust<-hclust(b.bray) #cluster
b.grp<- as.data.frame(cutree(b.clust, 3)) #define groups 
b.nmds<-metaMDS(diat[,6:62],distance="bray") 
b1 <- b.nmds$points[,1]
b2 <- b.nmds$points[,2]
bc<- as.data.frame(cbind(b1, b2, b.grp))
colnames(bc)[3]<- 'grp'
grp<- as.factor(bc$grp)

b.nmds.plot<- ggplot()
b.nmds.plot<- b.nmds.plot + geom_vline(x=0,colour="grey50")
b.nmds.plot<- b.nmds.plot + geom_hline(y=0,colour="grey50")
b.nmds.plot<- b.nmds.plot + geom_text(data = bc, aes(x=b1, y=b2, label = diat$Estimated_year, colour = grp), size = 7)
b.nmds.plot<- b.nmds.plot + theme_bw()
b.nmds.plot<- b.nmds.plot + labs(x= "NMDS 1", y= "NMDS 2")
b.nmds.plot<- b.nmds.plot + theme(axis.text.x = element_text(colour="black", size=16))
b.nmds.plot<- b.nmds.plot + theme(axis.text.y = element_text(colour="black", size=16))
b.nmds.plot<- b.nmds.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
b.nmds.plot<- b.nmds.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))
b.nmds.plot<- b.nmds.plot + ggtitle("NMDS - Counts")


################TEMPORAL BETA DIVERSITY USING VALVE DIATOM DATA (LAPERRIERE)

#Diatom valve data (counts)
diat<- read.csv("dauriat_subfossil_diatoms1998.csv") 

#Make each row a dataframe - begin at 1937.25

#T0 - DAR_Lap (oldest) - 1900
intT0<- as.data.frame(subset(diat, Interval_median == "37.25", drop=T))

#T1 - DAR_Lap  - 1900
intT1<- as.data.frame(subset(diat, Interval_median == "36.25", drop=T))

#T2 - DAR_Lap - 1920
intT2<- as.data.frame(subset(diat, Interval_median == "35.25", drop=T))

#T3 - DAR_Lap - 1920
intT3<- as.data.frame(subset(diat, Interval_median == "34.25", drop=T))

#T4 - DAR_Lap - 1930
intT4<- as.data.frame(subset(diat, Interval_median == "33.25", drop=T))

#T5 - DAR_Lap - 1930
intT5<- as.data.frame(subset(diat, Interval_median == "32.25", drop=T))

#T6 - DAR_Lap - 1939
intT6<- as.data.frame(subset(diat, Interval_median == "31.25", drop=T))

#T7 - DAR_Lap - 1939
intT7<- as.data.frame(subset(diat, Interval_median == "30.25", drop=T))

#T8 - DAR_Lap - 1950
intT8<- as.data.frame(subset(diat, Interval_median == "29.25", drop=T))

#T9 - DAR_Lap - 1950
intT9<- as.data.frame(subset(diat, Interval_median == "28.25", drop=T))

#T10 - DAR_Lap - 1960
intT10<- as.data.frame(subset(diat, Interval_median == "27.25", drop=T))

#T11 - DAR_Lap - 1960
intT11<- as.data.frame(subset(diat, Interval_median == "26.25", drop=T))

#T12 - DAR_Lap - 1965
intT12<- as.data.frame(subset(diat, Interval_median == "25.25", drop=T))

#T13 - DAR_Lap - 1965
intT13<- as.data.frame(subset(diat, Interval_median == "24.25", drop=T))

#T14 - DAR_Lap - 1970
intT14<- as.data.frame(subset(diat, Interval_median == "23.25", drop=T))

#T15 - DAR_Lap - 1970
intT15<- as.data.frame(subset(diat, Interval_median == "22.25", drop=T))

#T16 - DAR_Lap - 1974
intT16<- as.data.frame(subset(diat, Interval_median == "21.25", drop=T))

#T17 - DAR_Lap - 1974
intT17<- as.data.frame(subset(diat, Interval_median == "20.25", drop=T))

#T18 - DAR_Lap - 1977
intT18<- as.data.frame(subset(diat, Interval_median == "19.25", drop=T))

#T19 - DAR_Lap - 1977
intT19<- as.data.frame(subset(diat, Interval_median == "18.25", drop=T))

#T20 - DAR_Lap - 1980
intT20<- as.data.frame(subset(diat, Interval_median == "17.25", drop=T))

#T21 - DAR_Lap - 1980
intT21<- as.data.frame(subset(diat, Interval_median == "16.25", drop=T))

#T22 - DAR_Lap - 1982.5
intT22<- as.data.frame(subset(diat, Interval_median == "15.25", drop=T))

#T23 - DAR_Lap - 1982.5
intT23<- as.data.frame(subset(diat, Interval_median == "14.25", drop=T))

#T24 - DAR_Lap - 1985
intT24<- as.data.frame(subset(diat, Interval_median == "13.25", drop=T))

#T25 - DAR_Lap - 1985
intT25<- as.data.frame(subset(diat, Interval_median == "12.25", drop=T))

#T26 - DAR_Lap - 1987.5
intT26<- as.data.frame(subset(diat, Interval_median == "11.25", drop=T))

#T27 - DAR_Lap - 1987.5
intT27<- as.data.frame(subset(diat, Interval_median == "10.25", drop=T))

#T28 - DAR_Lap - 1990
intT28<- as.data.frame(subset(diat, Interval_median == "9.25", drop=T))

#T29 - DAR_Lap - 1990
intT29<- as.data.frame(subset(diat, Interval_median == "8.25", drop=T))

#T30 - DAR_Lap - 1992.5
intT30<- as.data.frame(subset(diat, Interval_median == "7.25", drop=T))

#T31 - DAR_Lap - 1992.5
intT31<- as.data.frame(subset(diat, Interval_median == "6.25", drop=T))

#T32 - DAR_Lap - 1995
intT32<- as.data.frame(subset(diat, Interval_median == "5.25", drop=T))

#T33 - DAR_Lap - 1995
intT33<- as.data.frame(subset(diat, Interval_median == "4.25", drop=T))

#T34 - DAR_Lap - 1997
intT34<- as.data.frame(subset(diat, Interval_median == "3.25", drop=T))

#T35 - DAR_Lap - 1997
intT35<- as.data.frame(subset(diat, Interval_median == "2.25", drop=T))

#T36 - DAR_Lap - 1999
intT36<- as.data.frame(subset(diat, Interval_median == "1.25", drop=T))

#T37 - DAR_Lap - 1999
intT37<- as.data.frame(subset(diat, Interval_median == "0.25", drop=T))


#Comparisons
#T0-T1
T0T1.int<- decompose.D2(intT0[,6:62], intT1[,6:62], den.type=2)  
T0T1.int.mat2<- as.data.frame(T0T1.int$mat2) 

#T1-T2
T1T2.int<- decompose.D2(intT1[,6:62], intT2[,6:62], den.type=2)  
T1T2.int.mat2<- as.data.frame(T1T2.int$mat2) 

#T2-T3
T2T3.int<- decompose.D2(intT2[,6:62], intT3[,6:62], den.type=2)  
T2T3.int.mat2<- as.data.frame(T2T3.int$mat2) 

#T3-T4
T3T4.int<- decompose.D2(intT3[,6:62], intT4[,6:62], den.type=2)  
T3T4.int.mat2<- as.data.frame(T3T4.int$mat2) 

#T4-T5
T4T5.int<- decompose.D2(intT4[,6:62], intT5[,6:62], den.type=2)  
T4T5.int.mat2<- as.data.frame(T4T5.int$mat2) 

#T5-T6
T5T6.int<- decompose.D2(intT5[,6:62], intT6[,6:62], den.type=2)  
T5T6.int.mat2<- as.data.frame(T5T6.int$mat2) 

#T6-T7
T6T7.int<- decompose.D2(intT6[,6:62], intT7[,6:62], den.type=2)  
T6T7.int.mat2<- as.data.frame(T6T7.int$mat2) 

#T7-T8
T7T8.int<- decompose.D2(intT7[,6:62], intT8[,6:62], den.type=2)  
T7T8.int.mat2<- as.data.frame(T7T8.int$mat2) 

#T8-T9
T8T9.int<- decompose.D2(intT8[,6:62], intT9[,6:62], den.type=2)  
T8T9.int.mat2<- as.data.frame(T8T9.int$mat2) 

#T9-T10
T9T10.int<- decompose.D2(intT9[,6:62], intT10[,6:62], den.type=2)  
T9T10.int.mat2<- as.data.frame(T9T10.int$mat2) 

#T10-T11
T10T11.int<- decompose.D2(intT10[,6:62], intT11[,6:62], den.type=2)  
T10T11.int.mat2<- as.data.frame(T10T11.int$mat2) 

#T11-T12
T11T12.int<- decompose.D2(intT11[,6:62], intT12[,6:62], den.type=2)  
T11T12.int.mat2<- as.data.frame(T11T12.int$mat2) 

#T12-T13
T12T13.int<- decompose.D2(intT12[,6:62], intT13[,6:62], den.type=2)  
T12T13.int.mat2<- as.data.frame(T12T13.int$mat2) 

#T13-T14
T13T14.int<- decompose.D2(intT13[,6:62], intT14[,6:62], den.type=2)  
T13T14.int.mat2<- as.data.frame(T13T14.int$mat2) 

#T14-T15
T14T15.int<- decompose.D2(intT14[,6:62], intT15[,6:62], den.type=2)  
T14T15.int.mat2<- as.data.frame(T14T15.int$mat2) 

#T15-T16
T15T16.int<- decompose.D2(intT15[,6:62], intT16[,6:62], den.type=2)  
T15T16.int.mat2<- as.data.frame(T15T16.int$mat2) 

#T16-T17
T16T17.int<- decompose.D2(intT16[,6:62], intT17[,6:62], den.type=2)  
T16T17.int.mat2<- as.data.frame(T16T17.int$mat2) 

#T17-T18
T17T18.int<- decompose.D2(intT17[,6:62], intT18[,6:62], den.type=2)  
T17T18.int.mat2<- as.data.frame(T17T18.int$mat2) 

#T18-T19
T18T19.int<- decompose.D2(intT18[,6:62], intT19[,6:62], den.type=2)  
T18T19.int.mat2<- as.data.frame(T18T19.int$mat2) 

#T19-T20
T19T20.int<- decompose.D2(intT19[,6:62], intT20[,6:62], den.type=2)  
T19T20.int.mat2<- as.data.frame(T19T20.int$mat2) 

#T20-T21
T20T21.int<- decompose.D2(intT20[,6:62], intT21[,6:62], den.type=2)  
T20T21.int.mat2<- as.data.frame(T20T21.int$mat2) 

#T21-T22
T21T22.int<- decompose.D2(intT21[,6:62], intT22[,6:62], den.type=2)  
T21T22.int.mat2<- as.data.frame(T21T22.int$mat2) 

#T22-T23
T22T23.int<- decompose.D2(intT22[,6:62], intT23[,6:62], den.type=2)  
T22T23.int.mat2<- as.data.frame(T22T23.int$mat2) 

#T23-T24
T23T24.int<- decompose.D2(intT23[,6:62], intT24[,6:62], den.type=2)  
T23T24.int.mat2<- as.data.frame(T23T24.int$mat2) 

#T24-T25
T24T25.int<- decompose.D2(intT24[,6:62], intT25[,6:62], den.type=2)  
T24T25.int.mat2<- as.data.frame(T24T25.int$mat2) 

#T25-T26
T25T26.int<- decompose.D2(intT25[,6:62], intT26[,6:62], den.type=2)  
T25T26.int.mat2<- as.data.frame(T25T26.int$mat2) 

#T26-T27
T26T27.int<- decompose.D2(intT26[,6:62], intT27[,6:62], den.type=2)  
T26T27.int.mat2<- as.data.frame(T26T27.int$mat2) 

#T27-T28
T27T28.int<- decompose.D2(intT27[,6:62], intT28[,6:62], den.type=2)  
T27T28.int.mat2<- as.data.frame(T27T28.int$mat2) 

#T28-T29
T28T29.int<- decompose.D2(intT28[,6:62], intT29[,6:62], den.type=2)  
T28T29.int.mat2<- as.data.frame(T28T29.int$mat2) 

#T29-T30
T29T30.int<- decompose.D2(intT29[,6:62], intT30[,6:62], den.type=2)  
T29T30.int.mat2<- as.data.frame(T29T30.int$mat2) 

#T30-T31
T30T31.int<- decompose.D2(intT30[,6:62], intT31[,6:62], den.type=2)  
T30T31.int.mat2<- as.data.frame(T30T31.int$mat2) 

#T31-T32
T31T32.int<- decompose.D2(intT31[,6:62], intT32[,6:62], den.type=2)  
T31T32.int.mat2<- as.data.frame(T31T32.int$mat2) 

#T32-T33
T32T33.int<- decompose.D2(intT32[,6:62], intT33[,6:62], den.type=2)  
T32T33.int.mat2<- as.data.frame(T32T33.int$mat2) 

#T33-T34
T33T34.int<- decompose.D2(intT33[,6:62], intT34[,6:62], den.type=2)  
T33T34.int.mat2<- as.data.frame(T33T34.int$mat2)

#T34-T35
T34T35.int<- decompose.D2(intT34[,6:62], intT35[,6:62], den.type=2)  
T34T35.int.mat2<- as.data.frame(T34T35.int$mat2)

#T35-T36
T35T36.int<- decompose.D2(intT35[,6:62], intT36[,6:62], den.type=2)  
T35T36.int.mat2<- as.data.frame(T35T36.int$mat2)

#T36-T37
T36T37.int<- decompose.D2(intT36[,6:62], intT37[,6:62], den.type=2)  
T36T37.int.mat2<- as.data.frame(T36T37.int$mat2)













#################PLOTTING FULL METAL PROFILES#################
##Plotting full metal profiles for Dauriat and Knob
#then univariate regression trees with y = metal EF and x = year (time)

metal<- read.csv("euk_metalEF_full.csv") #includes EFs for a suite of metals as opposed to just the composite. 


#Plot- works if use just total metal_EF
darmet.plot<- ggplot(metal, aes(x=Est_year, y=Total_EF, colour = Lake)) + geom_path() + coord_flip() + facet_wrap(~Lake)
darmet.plot<- darmet.plot + geom_point(aes(x=Est_year, y=Total_EF, shape=factor(Year_type)), size=4)
darmet.plot<- darmet.plot + scale_colour_manual(values= c("Dauriat " = "darkblue", "Knob" = "firebrick3"))  
darmet.plot<- darmet.plot + scale_shape_manual(values = c(16,1))
darmet.plot<- darmet.plot + labs(x = "Ann�e", y = "Facteur d'enrichissement m�tallique") + theme_bw() + theme(legend.position="none")
darmet.plot<- darmet.plot + theme(axis.text.x = element_text(colour="black", size=16))
darmet.plot<- darmet.plot + theme(axis.text.y = element_text(colour="black", size=16))
darmet.plot<- darmet.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
darmet.plot<- darmet.plot + theme(axis.title.y = element_text(size = rel(2), angle=90)) 

#PC plots
darmet2.plot<- ggplot(metal, aes(x=Est_year, y=Metal_PC1, colour = Lake)) + geom_path() + coord_flip() + facet_wrap(~Lake)
darmet2.plot<- darmet2.plot + geom_point(aes(x=Est_year, y=Metal_PC1, shape=factor(Year_type)), size=4)
darmet2.plot<- darmet2.plot + scale_colour_manual(values= c("Dauriat " = "darkblue", "Knob" = "firebrick3"))  
darmet2.plot<- darmet2.plot + scale_shape_manual(values = c(16,1))
darmet2.plot<- darmet2.plot + labs(x = "Ann�e", y = "Composante principale 1") + theme_bw() + theme(legend.position="none")
darmet2.plot<- darmet2.plot + theme(axis.text.x = element_text(colour="black", size=16))
darmet2.plot<- darmet2.plot + theme(axis.text.y = element_text(colour="black", size=16))
darmet2.plot<- darmet2.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
darmet2.plot<- darmet2.plot + theme(axis.title.y = element_text(size = rel(2), angle=90)) 



#Univariate regression
#urt<- rpart(Metal_EF~Est_year, method="anova", data=metal) 
#printcp(urt) 
#plot(urt)
#summary(urt)
#urt.prune<-  prune(urt, cp=urt$cptable[which.min(urt$cptable[,"xerror"]),"CP"])
#summary(urt.prune)

#par(mfrow=c(1,2)) 
#rsq.rpart(urt.prune)  

# R2      
#1-urt.prune$cptable[2,2]   
#plot(as.party(urt.prune), tp_args = list(id = FALSE))


#Hclust - metal_EFS
metal.bray<-vegdist(metal[,7:16], method="bray")

met.clust<-hclust(metal.bray) 
plot(met.clust, hang=-1)

met.clust2<-reorder.hclust(met.clust, metal.bray)
plot(met.clust2, hang=-1)

#Using ggdendro
ggdendrogram(met.clust2, rotate = FALSE, size = 2)

met.clust3 <- as.dendrogram(met.clust2)

clust.dat <- dendro_data(met.clust3, type = "rectangle")

clust.plot <- ggplot() + geom_segment(data = clust.dat$segments, aes(x = x, y = y, xend = xend, yend = yend)) 
clust.plot<- clust.plot + theme_bw()
clust.plot<- clust.plot + geom_text(data = clust.dat$labels, aes(x = x, y = y, label = metal$Est_year, colour = metal$Lake), size = 3, vjust = 2)
#clust.plot<- clust.plot + scale_colour_manual(values = c("DAR" = "darkblue", "KB" = "red2"), "Lake")


metal.bray<-vegdist(metal[,7:16], method="bray")
met.clust<-hclust(metal.bray)
met.grp<- cutree(met.clust, 3) #cutting into 3 groups 

met.nmds<-metaMDS(metal[,7:16],distance="bray") 

col<- c("grey", "red", "lightblue")
col[met.grp]

plot(met.nmds, type = "n", display = "sites")
points(met.nmds, col = col[met.grp], bg = col[met.grp], pch = 21)
legend("topright", legend = paste("Cluster", 1:3),
       col = col, pt.bg = col, bty = "n", pch = 21)

plot(met.nmds, type = "n", display = "sites")
text(met.nmds, col = col[met.grp], bg = col[met.grp], label=metal$Est_year)
legend("topright", legend = paste("Cluster", 1:3),
       col = col, pt.bg = col, bty = "n", pch = 21)

#Try with ggplot so can better visualize 
metal.bray<-vegdist(metal[,7:16], method="bray") #distance matrix 
met.clust<-hclust(metal.bray) #cluster
met.grp<- as.data.frame(cutree(met.clust, 3)) #define groups 
met.nmds<-metaMDS(metal[,7:16],distance="bray") 
met1 <- met.nmds$points[,1]
met2 <- met.nmds$points[,2]
metc<- as.data.frame(cbind(met1, met2, met.grp))
colnames(metc)[3]<- 'grp'
grp<- as.factor(metc$grp)

met.nmds.plot<- ggplot()
met.nmds.plot<- met.nmds.plot + geom_vline(x=0,colour="grey50")
met.nmds.plot<- met.nmds.plot + geom_hline(y=0,colour="grey50")
met.nmds.plot<- met.nmds.plot +geom_point(data = metc, aes(x=met1, y=met2, shape = metal$Lake, colour = grp), size = 7)
#met.nmds.plot<- met.nmds.plot + geom_text(data = metc, aes(x=met1, y=met2, label = metal$Est_year, colour = grp), size = 7)
met.nmds.plot<- met.nmds.plot + theme_bw()
met.nmds.plot<- met.nmds.plot + labs(x= "NMDS 1", y= "NMDS 2")
met.nmds.plot<- met.nmds.plot + theme(axis.text.x = element_text(colour="black", size=16))
met.nmds.plot<- met.nmds.plot + theme(axis.text.y = element_text(colour="black", size=16))
met.nmds.plot<- met.nmds.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
met.nmds.plot<- met.nmds.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))


met.nmds.plot<- ggplot()
met.nmds.plot<- met.nmds.plot + geom_vline(x=0,colour="grey50")
met.nmds.plot<- met.nmds.plot + geom_hline(y=0,colour="grey50")
met.nmds.plot<- met.nmds.plot +geom_point(data = metc, aes(x=met1, y=met2, shape = metal$Lake, colour = grp), size = 7)
met.nmds.plot<- met.nmds.plot + geom_text(data = metc, aes(x=met1, y=met2, label = metal$Est_year, colour = grp), size = 7)
met.nmds.plot<- met.nmds.plot + theme_bw()
met.nmds.plot<- met.nmds.plot + labs(x= "NMDS 1", y= "NMDS 2")
met.nmds.plot<- met.nmds.plot + theme(axis.text.x = element_text(colour="black", size=16))
met.nmds.plot<- met.nmds.plot + theme(axis.text.y = element_text(colour="black", size=16))
met.nmds.plot<- met.nmds.plot + theme(axis.title.x = element_text(size = rel(2), angle=00))
met.nmds.plot<- met.nmds.plot + theme(axis.title.y = element_text(size = rel(2), angle=90))

#DAR only
metal.DAR<- as.data.frame(subset(metal, Lake == "Dauriat ", drop=T))

DAR.bray<-vegdist(metal.DAR[,7:16], method="bray")

met.clust<-hclust(DAR.bray) 
plot(met.clust, hang=-1)

met.clust2<-reorder.hclust(met.clust, metal.bray)
plot(met.clust2, hang=-1)

DAR.bray<-vegdist(metal.DAR[,7:16], method="bray")
met.clust<-hclust(DAR.bray)
met.grp<- cutree(met.clust, 3) #cutting into 3 groups 

DAR.nmds<-metaMDS(metal.DAR[,7:16],distance="bray") 

col<- c("grey", "red", "lightblue")
col[met.grp]

plot(DAR.nmds, type = "n", display = "sites")
points(DAR.nmds, col = col[met.grp], bg = col[met.grp], pch = 21)
legend("topright", legend = paste("Cluster", 1:3),
       col = col, pt.bg = col, bty = "n", pch = 21)

plot(DAR.nmds, type = "n", display = "sites")
text(DAR.nmds, col = col[met.grp], bg = col[met.grp], label=metal.DAR$Est_year)
legend("topright", legend = paste("Cluster", 1:3),
       col = col, pt.bg = col, bty = "n", pch = 21)


#KB only
metal.KB<- as.data.frame(subset(metal, Lake == "Knob", drop=T))

KB.bray<-vegdist(metal.KB[,7:16], method="bray")

met.clust<-hclust(KB.bray) 
plot(met.clust, hang=-1)

met.clust2<-reorder.hclust(met.clust, KB.bray)
plot(met.clust2, hang=-1)

KB.bray<-vegdist(metal.KB[,7:16], method="bray")
met.clust<-hclust(KB.bray)
met.grp<- cutree(met.clust, 2) #cutting into 2 groups 

KB.nmds<-metaMDS(metal.KB[,7:16],distance="bray") 

col<- c("grey", "red", "lightblue")
col[met.grp]

plot(KB.nmds, type = "n", display = "sites")
points(KB.nmds, col = col[met.grp], bg = col[met.grp], pch = 21)
legend("topright", legend = paste("Cluster", 1:3),
       col = col, pt.bg = col, bty = "n", pch = 21)

plot(KB.nmds, type = "n", display = "sites")
text(KB.nmds, col = col[met.grp], bg = col[met.grp], label=metal.KB$Est_year)
legend("topright", legend = paste("Cluster", 1:3),
       col = col, pt.bg = col, bty = "n", pch = 21)



######################################################
##Come back to everything below here###############

#######################NMDS##########################

#####################################################

################COMMON RICHNESS INTER SAMPLES (OTUs, Family)#########
a <-read.table("Nyme_cleaned_EMBL_nosingA+9r_48e5000r_families.txt",h=T)
data.frame(colnames(a))

b <- shary_nbtaxa(data=a,nb=48)
write.table(b,file="Nyme_cleaned_EMBL_nosingA+9r_48e5000r_nbcfam.txt",sep="\t")
b <- shary_perctaxa(data=a,nb=48)
write.table(b,file="Nyme_cleaned_EMBL_nosingA+9r_48e5000r_percfam.txt",sep="\t")


##############COMMON RICHNESS PAIRED-SAMPLES###################
a <-read.table("Nyme_cleaned_EMBL_nosingA_48e5497r_EMBL3_nbreads.txt",h=T)
data.frame(colnames(a))

z <-list()
for(i in 2:49)
{z[[i]]<-as.numeric(a[,i]>0)}
a2 <-data.frame(matrix(unlist(z),ncol=48))

z2 <-list()
for(i in 1:24)
{z2[[i]]<-cbind(data.frame(rowSums(a2[,(2*i-1):(2*i)])))}
a3 <-data.frame(z2)

z3 <- list()
for (i in 1:24)
{z3[[i]] <- ifelse(a3[,i] <= 1, 0, 1)}
a4 <- data.frame(matrix(unlist(z3),ncol=24))

a5 <- cbind(a[,1], a4, a[,50:length(a)])
colnames(a5) <- c("Family","N13_7","N13_10","N13_13","N13_16","N13_19","N13_22","N13_25","N13_28","N13_31","N13_34","N13_37","N13_40","N7_1","N7_4","N7_7","N7_10","N7_13","N7_16","N7_19","N7_22","N7_25","N7_28","N7_31","N7_34", colnames(a[,50:length(a)]))
data.frame(colnames(a5))

b <- taxy_nbreads(data=a5,number_of_samples=24,column_of_taxo_level=28,name_of_taxo_level="EMBL3")

write.table(b,file="Nyme_cleaned_EMBL_nosingA_48e5497r_nbfam_com.txt",sep="\t",row.names=F)

a<-read.table("Nyme_nosingB_common_fam_maingroups.txt",h=T)
data.frame(colnames(a))
barplot(as.matrix(t(a[,7:10])),cex.names=1,names.arg=toupper(a[,1]),  col=c("orangered", "gray","royalblue","gray"),las=2, horiz=T,width=c(5))

######################DOMINANCE RARETE######################
a <- read.table("Nyme_cleaned_EMBL_nosingA+1r_48e5497r.txt",h=T)
data.frame(colnames(a))

z <- list()
for (i in 2:49)
{z[[i]] <- a[,i]*100/sum(a[,i])}
a2 <- data.frame(matrix(unlist(z), ncol=48))
a2
a3 <- cbind(a[,1],a2, a[,50:length(a)])
write.table(a3,file="Nyme_cleaned_EMBL_nosingA+1r_48e5497r_perc.txt",sep="\t",row.names=F)


## number of taxa
z2 <- list()
for(i in 1:48)
{z2[[i]] <-  nrow(subset(a2,a2[,i] == 0))}

b0 <- unlist(z2)
b0.1 <- unlist(z2)
b1 <- unlist(z2)
b10 <- unlist(z2)
bsup10 <- unlist(z2)
bsup20 <- unlist(z2)
bsup30 <- unlist(z2)
bsup40 <- unlist(z2)
bsup50 <- unlist(z2)

A <- b0
B <- b0.1 -b0
C <- b1 - b0.1
D <- b10 - b1
E <- bsup10 - bsup20
F <- bsup20 - bsup30
G <- bsup30 - bsup40
H <- bsup40 - bsup50
I <- bsup50
AA <- cbind(A,B,C,D,E,F,G,H,I)
write.table(AA,file="Nyme_cleaned_EMBL_nosingA+1r_48e5497r_domrare.txt",sep="\t",row.names=F)

a <- read.table("nyme_domrare.txt",h=T)
data.frame(colnames(a))

colA <- c("gray", "orangered","green4","indianred")
barplot(as.matrix(a[1:8,2:49]), beside=F, col=colA)
        ,col=a1,las=2)
   
   
a <- read.table("nyme_SIMPER_EMBL3.txt")
data.frame(colnames(a))

colA <- c(rep("green4",2), rep("royalblue", 4), rep("chocolate4",2), rep("orangered",6))
CI = c(1,0,1,0,1,0,1,0,1,0,1,0,1,0)
par(mfcol=c(1,12))
barplot(as.matrix(t(a[,2])), beside=T, col=colA, horiz=T, space=c(CI), xlim=c(0,5000))
barplot(as.matrix(t(a[,3])), beside=T, col=colA, horiz=T, space=c(CI), xlim=c(0,5000))
barplot(as.matrix(t(a[,4])), beside=T, col=colA, horiz=T, space=c(CI), xlim=c(0,5000))
barplot(as.matrix(t(a[,5])), beside=T, col=colA, horiz=T, space=c(CI), xlim=c(0,5000))
barplot(as.matrix(t(a[,6])), beside=T, col=colA, horiz=T, space=c(CI), xlim=c(0,5000))
barplot(as.matrix(t(a[,7])), beside=T, col=colA, horiz=T, space=c(CI), xlim=c(0,5000))
barplot(as.matrix(t(a[,8])), beside=T, col=colA, horiz=T, space=c(CI), xlim=c(0,5000))
barplot(as.matrix(t(a[,9])), beside=T, col=colA, horiz=T, space=c(CI), xlim=c(0,5000))
barplot(as.matrix(t(a[,10])), beside=T, col=colA, horiz=T, space=c(CI), xlim=c(0,5000))
barplot(as.matrix(t(a[,11])), beside=T, col=colA, horiz=T, space=c(CI), xlim=c(0,5000))
barplot(as.matrix(t(a[,12])), beside=T, col=colA, horiz=T, space=c(CI), xlim=c(0,5000))
barplot(as.matrix(t(a[,13])), beside=T, col=colA, horiz=T, space=c(CI), xlim=c(0,5000))

#########################################