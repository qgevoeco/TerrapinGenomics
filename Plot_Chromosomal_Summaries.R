#!/usr/bin/env Rscript

#This script takes a single chromosome VCF as input and generates some windowed summary statistics using the popgenome package
#Include three command line arguments
# 1) path to a folder containing the desired VCF file
# 2) a text file with sample IDs in first column with title "ind", Poppulation designations in the second column with title "pop"
# 3) path to output directory for tables and figures

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
  stop("Supply four arguments in the following order. 1) path to a folder containing the desired VCF file. 2) a text file with sample IDs in first column with title [ind], Poppulation designations in the second column with title [pop] 3) path to output directory for tables and figures. 4) String indicating which Chromosome is being analyzed (e.g. Chr5)", call.=FALSE)
}

#################################
#VCF statistics
#################################
# Load packages
library(tidyverse)
library(PopGenome)
library(ggplot2)

############################################################
#Read in data and add population designations
############################################################
#Read in our genetic data, get the number of sites on the chromosome
GENOME_OBJ <- readData(args[1], format = "VCF", include.unknown = TRUE, FAST = TRUE) #Creates a GENOME object
sum<-get.sum.data(GENOME_OBJ)
sum<-as.data.frame(sum)
N_Sites<-sum$n.sites

#Read in the population data
Info<- read_delim(args[2], delim = "\t")
populations <- split(Info$ind, Info$pop) #Split the data into lists containing the indivs in each population
Genome_with_pops <- set.populations(GENOME_OBJ, populations, diploid = T) #Add the populations to the GENOME object
#Verify that the populations look correct
Genome_with_pops@populations

####################################
#Set up sliding window analyses
####################################
#Set chromosome size
Chrom_Size <- N_Sites

#Set window size and window jump- As a default we'll go with 100k windows that jump 25k, so these windows will be overlapping
window_size <- 100000
window_jump <- 25000

# use seq to find the start points of each window
window_start <- seq(from = 1, to = Chrom_Size, by = window_jump)
# add the size of the window to each start point 
window_stop <- window_start + window_size
# no windows start before the end of chromosome 8
sum(window_start > Chrom_Size)
# but some window stop positions do occur past the final point
sum(window_stop > Chrom_Size)
# remove windows from the start and stop vectors
window_start <- window_start[which(window_stop < Chrom_Size)]
window_stop <- window_stop[which(window_stop < Chrom_Size)]
# save as a data.frame
windows <- data.frame(start = window_start, stop = window_stop, 
                      mid = window_start + (window_stop-window_start)/2)

##########################################################################################################
#Conduct sliding window analyses to get estimates of Pi within pops and Dxy + Fst among populations 
##########################################################################################################
#Make a sliding window dataset (make 100% sure this is the same as above)
Separate_Pop_sw <- sliding.window.transform(Genome_with_pops, width = 100000, jump = 25000, type = 2)
Separate_Pop_sw <- diversity.stats(Separate_Pop_sw, pi = TRUE) #Adds nucleotide diversity (pi + dxy) to the object
Separate_Pop_sw <- F_ST.stats(Separate_Pop_sw, mode = "nucleotide") #Adds Fst to the object

#Extract stats for visualization into the nd dataframe for dxy and pi, fst frame for Fst
nd <- Separate_Pop_sw@nuc.diversity.within/100000 #Divides the number of differences by the width of the window to get pi

# make population name vector that we'll use later (alphabetical order please)
pops <- c("CedarPoint", "Dauphin", "Mississippi")

# set population names
colnames(nd) <- paste0(pops, "_pi")

#Extract fst values
fst <- t(Separate_Pop_sw@nuc.F_ST.pairwise) #transpose the paiwise Fst table from the sliding window object
fst[is.na(fst)] <- 0 #Set NA values of Fst to zero
fst[fst < 0] <- 0 #Set negative Fst values to zero

#Extract Dxy values
dxy <- get.diversity(Separate_Pop_sw, between = T)[[2]]/100000

# get column names 
x <- colnames(fst)

# swap in our specific population names for the generic pop1,2,3 stuff
x <- sub("pop1", pops[1], x)
x <- sub("pop2", pops[2], x)
x <- sub("pop3", pops[3], x)

# replace forward slash with an underscore
x <- sub("/", "_", x)

#Change column names
colnames(fst) <- paste0(x, "_fst")
colnames(dxy) <- paste0(x, "_dxy")

#Combine the fst, pi, and dxy dataframes with the data frame containing position along the chromosome
All_data <- data.frame(windows, nd, fst, dxy)
CP_pi<-mean(All_data$CedarPoint_pi)
DI_pi<-mean(All_data$Dauphin_pi)
MS_pi<-mean(All_data$Mississippi_pi)

# select nucleotide diversity data and calculate means
All_data %>% select(contains("_pi")) %>% summarise_all(mean)

# gather the data
pi_g <- All_data %>% select(contains("_pi")) %>% gather(key = "species", value = "pi")

##################################################################
# Make some pretty plots along the chromosome
##################################################################
# select data of interest
hs <- All_data %>% select(mid, CedarPoint_pi, Dauphin_pi, CedarPoint_Dauphin_fst, CedarPoint_Dauphin_dxy)
hs_g <- gather(hs, -mid, key = "stat", value = "value")
hs_g <- hs_g[!duplicated(hs_g), ]
# construct a plot with facets
a <- ggplot(hs_g, aes(mid, value, colour = stat)) + geom_line()
a <- a + facet_grid(stat~., scales = "free_y")
a <- a + xlab("Position (Mb)")
a + theme_light() + theme(legend.position = "none")
#Now, write the plot to a .png file
setwd(args[3])
png("Summmary_Plot.png")
print(a)
dev.off()



