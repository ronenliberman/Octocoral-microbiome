## 16S Microbiome Analyses - Molecular Microbiology & Genomics Lab ##
## This R-Script includes routine microbiome analyses conducted using
## the 16S rRNA gene (V4 region (primers 515F/806R)).
## Initial data cleaning starts from QIIME2, so you will need the following 
## files from the QIIME2 resulting files in your desired working folder: 
## feature.csv and taxonomy.csv (both are originally tsv files but use Excel to convert to csv)
## Feel free to change anything in this script to fit your project analyses!
## Please make sure to keep this original copy and to save any edited versions
## with your project data!

## This is a link to the entire color chart used by R: https://r-charts.com/colors/
## You can also run the code: colors() and R will give you a list of the color names.
## I personally use the link since you can see the colors.

## Created by Paisley Samuel - Summer 2023 ##

###### Table of Contents/Order of Processes ######
### 1. List of packages that are used throughout the script 
### 2. Setting the working directory and seed 
### 3. Producing the abundance data using raw sequencing data
### 4. Checking the produced abundance data for the batch effect and correcting if necessary
### 5. Generating a rarefaction curve of the original/adjusted abundance data
### 6. Reproducing the abundance data using the batch-corrected sequencing data
### 7. Taxonomy analyses on abundance data
### 8. Creating a table of calculated alpha diversity metrics
### 9. Performing statistical analyses between metadata variables and alpha diversity metrics
### 10. Creating a Bray-Curtis distance matrix and performing beta diversity statistical analyses
### 11. Creating nMDS plots using BC distance matrix 
### 12. Performing CCA analysis and creating CCA plots using 

###### List of Packages Used ######
## Here is a list of all the packages used throughout this R script.
## Take the necessary steps to install any packages that you do not have
## installed already on your computer (may have to use install.packages()
## or Bioconductor depending on the package)

## install qiime2R
#if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
#devtools::install_github("jbisanz/qiime2R")

library(ggplot2)
library(MBiome16S) # created specifically for 16S V4 region analyses in Lopez lab
library(microbiome)
library(MMUPHin)
library(multcompView)
library(pgirmess)
library(phyloseq)
library(psych)
library(tidyverse)
library(vegan)
library(qiime2R) 

###### Set Working Directory and Seed - ALWAYS RUN THIS FIRST! ######
setwd_seed("~/Postdocing_NSU/JOE LOPEZ Lab/Lab WORK/octocoral 16S/Merged_June_October_24/Merged_analysis/Octocoral-microbiome/raw files/", 23)
###### Generating Rarefaction Curve #######

## assigning the batch-corrected or original feature table to its own variable; MAKE SURE TO READ ROWNAMES

rardat<-read.csv("~/Postdocing_NSU/JOE LOPEZ Lab/Lab WORK/octocoral 16S/Merged_June_October_24/Merged_analysis/Octocoral-microbiome/raw files/feature-table.csv" ,
header=TRUE, row.names=1, sep=',')

# samples are in columns and need to be in the rows so we need to flip or transpose the file
# transpose the data to rows; transposing changes the data frame to a matrix!
trans.rardat <- t(rardat)
# check file to make sure it worked
trans.rardat[1:5,1:5] #shows rows 1 through 5 and the samples should now be the rows
# assign the transformed data matrix into main data variable
rardat <- trans.rardat
# change back into data frame instead of matrix
rardat <-as.data.frame(rardat)
#check data file to make sure it looks okay
head(rardat)

rowSums(rardat) #sums the value of each row in the data frame;
# this shows the total sequencing reads for each sample

#looking for sampoles with very low count 
# Identify samples with very low reads (between 0 and 10)
total_reads <- rowSums(rardat)
total_reads[total_reads <= 10] 
#' now I can delete samples: 
#' EFLE.DC1.3.2024 EFLE.DC3.1.2024 EFLE.DC3.3.2024 
#' Now I created a new feawture table called "feature-table.UPD.csv" and i will be using it instead of the original one.

## Creating the rarefaction curve
# count the number of species within each sample
S <- specnumber(rardat)

raremax <- min(rowSums(rardat)) # takes the sample with the lowest number of sequencing reads

## Plotting the rarefaction curves
# ** auto removes samples that have no reads **

# creating color palette
col <- c("darkred", "forestgreen", "hotpink", "blue")  # feel free to edit the colors

# keep it between 3 and 5 colors
grp <- factor(sample(seq_along(col), nrow(rardat), replace = TRUE))
cols <- col[grp]

# Creating rarefaction curve
# create the curve to estimate where the inflection point (i.e. the point at which
# most of the lines being to plateau) lies, then assign that value to the
# variable "inflection" below
png("RarefactionCurve_Jun24_Oct24run.png", res = 300, width = 3900, height = 2700) #creates a file in your working folder to house the graph when it populates below
rarecurve(rardat, step = 500, col = cols, label = TRUE,
          main="Rarefraction of octocoral sequencing of 4 species sampled in summer 23 and 24", cex= 0.35, cex.axis= 0.95, cex.lab= 1, xlim=c(0,175000),
          xlab = "# of Sequencing Reads", ylab = "# of ASVs")
inflection <- 20000   # insert your estimated inflection point value
abline(v = inflection, col="black", lwd=1.4) # creates the inflection line at specified value
# feel free to change the line color and thickness (col and lwd)
graphics.off() #saves the graph to the file

# generating a taxonomy files based on phylogeny 

taxonomy.raw<-read_qza("taxonomy.qza") #using raw tax file from qiime2
head(taxonomy$data)

taxonomy<-parse_taxonomy(taxonomy.raw$data)
taxonomy$Kingdom<-gsub("d__","", taxonomy$Kingdom)

write.csv(file="new_taxonomy.csv", taxonomy)

###### Producing Relative Abundance Data ######

#abund_metadata <- calculate_abund.metadata("~/Postdocing_NSU/JOE LOPEZ Lab/Lab WORK/octocoral 16S/Merged_June_October_24/Merged_analysis/Octocoral-microbiome/raw files/feature-table.csv" ,
                    #"~/Postdocing_NSU/JOE LOPEZ Lab/Lab WORK/octocoral 16S/Merged_June_October_24/Merged_analysis/Octocoral-microbiome/raw files/metadata.csv")

abund_metadata <- calculate_abund.metadata("feature-table.UPD.csv" ,"metadata.csv")

# your workspace, simply assign it to its own variable
# Example saving dat and metadata
feature <- abund_metadata$dat
metadata <- abund_metadata$metadata
dat.01per <- abund_metadata$dat.01per
dat.ra <- abund_metadata$dat.ra

head(metadata)

#tax= read.csv("new_taxonomy.csv")

## Basic sequencing statistics (USING THE "physeq" OBJECT IN YOUR PHYLOSEQ OBJECT LIST VARIABLE)
# First, assign physeq and physeq_transform to their own variables
physeq_objs <- create_phyloseq(dat.01per, "new_taxonomy.csv", "metadata.csv")

# abundance = insert the abundance table you created (dat.01per)
# "metadata.csv" = insert your metadata .csv file

## Basic sequencing statistics (USING THE "physeq" OBJECT IN YOUR PHYLOSEQ OBJECT LIST VARIABLE)
# First, assign physeq and physeq_transform to their own variables
physeq <- physeq_objs$physeq
physeq_transform <- physeq_objs$physeq_transform

# Check number of sequencing reads observed in each sample
sample_sums(physeq)

# Calculate the total sequencing reads of your data
sum(sample_sums(physeq)) #11796103

# Calculate the average number sequencing reads
mean(sample_sums(physeq)) #77605.94

# Find the lowest number of sequencing reads in your data
min(sample_sums(physeq)) #8039

# Find the highest number of sequencing reads in your data
max(sample_sums(physeq)) #232408

# Calculate the standard deviation of sequencing reads between samples
sd(sample_sums(physeq)) #258994

# Find the total amount of ASVs within your data
ntaxa(physeq) #16334

## Retrieves the unique taxonomic ranks observed in the data set
## [#] = rank (starting from Domain and onward DPCOFGS)
get_taxa_unique(physeq, taxonomic.rank=rank_names(physeq)[], errorIfNULL=TRUE)
#Unique Domains [1] =
#Unique Phyla [2] =
#Unique Classes [3] =
#Unique Orders [4] =            #SKIPPED!
#Unique Families [5] =
#Unique Genera [6] =
#Unique Species [7] =

## Aggregating by Taxonomic level
# This function allows you to aggregate the taxonomy based on Taxonomic
# level, gives you the counts of each level, and saves as a CSV file

#glimpse(physeq_transform)
# phyloseq_object = insert the transformed phyloseq object you created (physeq_transform)
# "tax_level" = insert the taxonomic level you will like to aggregate on


# phyloseq_object = insert the transformed phyloseq object you created (physeq_transform)
# "tax_level" = insert the taxonomic level you will like to aggregate on
# n = insert the number of taxa you want (e.g. top 5, top 10, etc.)

# For example for top 5 Classes, the code would read:
# top.class.names <- top_taxa(physeq_transform,"Class",5)


# Next, subset your phyloseq object to only include the top taxa you just specified
# cuts down the phyloseq object to only the top n
# Replace the "taxonomic level" with the level you want WITHOUT QUOTES!


top.ten.class <- MBiome16S::top_taxa(physeq_transform, "Class", 10)
top.ten.class

top.ten.order <- MBiome16S::top_taxa(physeq_transform, "Order", 10)
top.ten.order

top.ten.family <- MBiome16S::top_taxa(physeq_transform, "Family", 10)
top.ten.family

top.ten.gen <- MBiome16S::top_taxa(physeq_transform, "Genus", 10)
top.ten.gen

##Exporting taxonomy 
## Creating stacked bar plots from your phyloseq object subset
# Change any taxonomic names to the level that you need ("fill=" and "colour=")
top_class<- subset_taxa(physeq_transform, Class %in% names(top.ten.class))
topTaxplot_class <- plot_bar(top_class, x="Sample", y="Abundance", fill="Class")
topTaxplot_class <- topTaxplot_class +
  geom_bar(aes(fill=Class, colour=Class), stat="identity", position="fill", width = 0.9) +   #width=0.96 removes any space between bars
  ggtitle("Top 10 Classes") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.28))+   #vjust= moves the x-axis text labels
  theme(plot.title = element_text(color="black", size=14, face="bold.italic", hjust = 0.5)) +   #hjust= 0.5 centers the title
  theme(legend.title = element_text(face="italic"))
topTaxplot_class 

### order 
top_order<- subset_taxa(physeq_transform, Order %in% names(top.ten.order))

## Creating stacked bar plots from your phyloseq object subset
# Change any taxonomic names to the level that you need ("fill=" and "colour=")

topTaxplot_order <- plot_bar(top_order, x="Sample", y="Abundance", fill="Order")
topTaxplot_order <- topTaxplot_order +
  geom_bar(aes(fill=Order, colour=Order), stat="identity", position="fill", width = 0.9) +   #width=0.96 removes any space between bars
  ggtitle("Top 10 Orders") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.28))+   #vjust= moves the x-axis text labels
  theme(plot.title = element_text(color="navyblue", size=14, face="bold.italic", hjust = 0.5)) +   #hjust= 0.5 centers the title
  theme(legend.title = element_text(face="italic"))
topTaxplot_order 

#Family

top_family<- subset_taxa(physeq_transform, Family %in% names(top.ten.family))
topTaxplot_family<- plot_bar(top_family, x="Sample", y="Abundance", fill="Family")
topTaxplot_family <- topTaxplot_family +
  geom_bar(aes(fill=Family, colour=Family), stat="identity", position="fill", width = 0.9) +   #width=0.96 removes any space between bars
  ggtitle("Top 10 Family") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.28))+   #vjust= moves the x-axis text labels
  theme(plot.title = element_text(color="black", size=14, face="bold.italic", hjust = 0.5)) +   #hjust= 0.5 centers the title
  theme(legend.title = element_text(face="italic"))
topTaxplot_family 


###### Alpha Diversity - Measures ######
## Alpha diversity: the species richness that occurs within a given area within a region
## that is smaller than the entire distribution of the species (Moore, 2013)

### I use the relative abundance or dat.01per data
adivmeasures(dat.01per)
## This function returns one dataframe that contains the alpha diversity 
## measure metrics for each sample


# Merge diversity with metadata table and export as csv
# DONE IN EXCEL BEFORE RUNNING NEXT CODE: rename the first column of both data
# tables as the SAME NAME. e.g. both should have the sample column labelled as "Sample"
diversitybysample <- read.csv("AlphaDiversity.csv", row.names = 1)

#met <- read.csv("metadata.csv", row.names = 1)

rownames(metadata)

#comparing row names of my alpha diversity data and metadata
# Extract row names
div_row <- rownames(diversitybysample)
metadata_row <- rownames(metadata)
# Find common row names
common_row_names <- intersect(div_row, metadata_row)

# Subset data frames to keep only common rows
diversitybysample <- diversitybysample[common_row_names, , drop = FALSE]
met <- metadata[common_row_names, , drop = FALSE]

adivmet <- cbind(diversitybysample,met)
head(adivmet)
write.csv(adivmet,"Metadata-Diversity.csv")

###### Alpha Diversity Statistics ######

# load in your metadata that has the alpha diversity indices included
Alpha_data <- read.csv("Metadata-Diversity.csv", header = TRUE, row.names = 1)
head(Alpha_data)

#### Testing Statistical Significance
## Normality - Shapiro Test (only done on NUMERIC data)
## p <= 0.05 = H0 REJECTED -> DATA IS NOT NORMAL
## p > 0.05 = H0 ACCEPTED -> DATA IS NORMAL

## If data is NOT normal the first time, try transforming the data using log
## and sqrt and retest for normality.

#Alpha Diversity Variables - Test for Normality
shapiro.test(Alpha_data$S) #not normal
shapiro.test(Alpha_data$N) #not normal
shapiro.test(Alpha_data$H) #not Normal
shapiro.test(Alpha_data$J) #not Normal
shapiro.test(Alpha_data$inv.D) #not normal

str(Alpha_data)

#format as factor for categorical factors
Alpha_data$SPP <- as.factor(Alpha_data$SPP)
Alpha_data$Site <- as.factor(Alpha_data$Site)
Alpha_data$Year <- as.factor(Alpha_data$Year)
Alpha_data$Region <- as.factor(Alpha_data$Region)

str(Alpha_data)
#Site

### CONTINUE HERE IF DATA IS NOT NORMAL AND TRANSFORMATIONS DID NOT WORK
# Kruskal Wallis: Nonparametric Data (not normal)
# Pairwise Wilcox Test - calculate pairwise comparisons between group levels
#     with corrections for multiple testing (non-parametric)

#Species richness
S_krusk_Species <- nonp_kruskal(Alpha_data$S, Alpha_data$SPP)  #significant
S_krusk_Year <- nonp_kruskal(Alpha_data$S, Alpha_data$Year) #significant , but very close to 0.05
S_krusk_Region <- nonp_kruskal(Alpha_data$S, Alpha_data$Region) #significantbut very close to 0.05
S_krusk_Site <- nonp_kruskal(Alpha_data$S, Alpha_data$Site) #significant, but most sites are similar except for DC1-BC5 which anyway dont have all species

#Shannon Diversity Index (H)
H_krusk_Species <- nonp_kruskal(Alpha_data$H, Alpha_data$SPP) #MMUR is different
H_krusk_Year <- nonp_kruskal(Alpha_data$H, Alpha_data$Year) #significant , but very close to 0.05
H_krusk_Region <- nonp_kruskal(Alpha_data$H, Alpha_data$Region) #significantbut very close to 0.05
H_krusk_Site <- nonp_kruskal(Alpha_data$H, Alpha_data$Site) #No difference

#Species Evenness (J)

J_krusk_Species <- nonp_kruskal(Alpha_data$J, Alpha_data$SPP) #MMUR is different
J_krusk_Year <- nonp_kruskal(Alpha_data$J, Alpha_data$Year) #significant , but very close to 0.05
J_krusk_Region <- nonp_kruskal(Alpha_data$J, Alpha_data$Region) #significantbut very close to 0.05
J_krusk_Site <- nonp_kruskal(Alpha_data$J, Alpha_data$Site) #No difference

#number of individuals
inv.D_krusk_Species <- nonp_kruskal(Alpha_data$inv.D, Alpha_data$SPP) #MMUR is different
inv.D_krusk_Site <- nonp_kruskal(Alpha_data$inv.D, Alpha_data$Site) #No difference

# Looking into the Kruskal-Wallis results for each diversity metric
# The important results to be noted are the KW p-values and the letter comparisons.
# You can also take a look at the individual comparison p-values by looking into the pairwise table.

S_krusk_Species$`Kruskal-Wallis Results` ## Kruskal-Wallis chi-squared = 44.745, df = 3, p-value = 1.048e-09
S_krusk_Species$pairwise
S_krusk_Species$`letter-comparisons` # Bast = SW, EFLE = MMUR, any other combination is not significant 

### You would then repeat these steps for the 4 other alpha diversity metrics!

### Plotting boxplots of alpha diversity by specified variable
## NOTES: You can replace "Year" with your specified variable
##        Adding text to your graph is OPTIONAL but if you are adding it then
##        you'll have to play around with their coordinates, what they say, size, and color.

# Creating the plots for different factors
# Load necessary libraries
library(ggplot2)
library(patchwork)

# Define a colorblind-friendly palette
cb_palette <- scale_fill_manual(values = c("#D55E00", "#009E73", "#F0E442", "#0072B2", "#CC79A7"))

# Adjusted Species Richness (S) plot
S_krusk_Species$`letter-comparisons`

p_S <- ggplot(data = Alpha_data, aes(x = SPP, y = S)) +
  geom_violin(aes(fill = SPP), alpha = 0.6) +
  annotate("text", x = 1, y = 2500, label = "a", size = 5) +
  annotate("text", x = 2, y = 2500, label = "b", size = 5) +
  annotate("text", x = 3, y = 2000, label = "c", size = 5) +
  annotate("text", x = 4, y = 700, label = "c", size = 5) +
  annotate("text", x = 5, y = 1100, label = "ac", size = 5) +
  labs(title = "Species Richness (S)", y = "Species Richness", x = NULL) +
  theme_bw() +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  cb_palette

# Adjusted Shannon Diversity Index (H) plot
p_H <- ggplot(data = Alpha_data, aes(x = SPP, y = H)) +
  geom_violin(aes(fill = SPP), alpha = 0.6) +
  annotate("text", x = 1, y = 6.7, label = "a", size = 5) +
  annotate("text", x = 2, y = 6.7, label = "b", size = 5) +
  annotate("text", x = 3, y = 6, label = "c", size = 5) +
  annotate("text", x = 4, y = 4, label = "c", size = 5) +
  annotate("text", x = 5, y = 6.3, label = "ab", size = 5) +
  labs(title = "Shannon Diversity Index (H)", y = "Shannon Diversity", x = NULL) +
  theme_bw() +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  cb_palette

# Adjusted Species Evenness (J) plot
p_J <- ggplot(data = Alpha_data, aes(x = SPP, y = J)) +
  geom_violin(aes(fill = SPP), alpha = 0.6) +
  annotate("text", x = 1, y = 0.95, label = "a", size = 5) +
  annotate("text", x = 2, y = 0.95, label = "b", size = 5) +
  annotate("text", x = 3, y = 0.85, label = "c", size = 5) +
  annotate("text", x = 4, y = 0.75, label = "c", size = 5) +
  annotate("text", x = 5, y = 0.9, label = "ab", size = 5) +
  labs(title = "Species Evenness (J)", y = "Species Evenness", x = "Species") +
  theme_bw() +
  theme(legend.position = "none") +
  cb_palette

# Combine the plots using patchwork
combined_plot <- (p_S / p_H / p_J) +
  plot_layout(guides = "collect") &
  theme(plot.title = element_text(size = 14, face = "bold"),
        plot.caption = element_text(size = 10))

# Display the combined plot
print(combined_plot)

ggsave(
  filename = "SPP_comparison_combined_Alpha diversity.jpeg", # Output file name
  plot = combined_plot,           # Plot object
  device = "jpeg",                # File format
  width = 10,                     # Width of the image in inches
  height = 12,                    # Height of the image in inches
  dpi = 300                       # Resolution in dots per inch
)


#################

bc.dist <- vegdist(dat.01per, method = "bray") #creating distance matrix 

dis.SPP <- betadisper(bc.dist, metadata$SPP) #calculates the dispersion (variances) within each group

permutest(dis.SPP, pairwise = T, permutations = 999) #determines if the variances differ by groups

# permutest significant - USING ANOSIM
## IF permutest NOT SIGNIFICANT = USE PERMANOVA (adonis2)

anosim(bc.dist,  metadata$SPP, permutations = 999)

###### Beta Diversity - nMDS plots ######
# Creating nMDS plot - 2D
nmds2d <- metaMDS(bc.dist,k=2,autotransform = F,trymax=50) # you will have to extract bc.dist from the previous analysis (object betadiv)
nmds2d 

# Dimensions: 2 
# Stress: 0.2071628 -> your stress it a bit high

stressplot(nmds2d) 

#Shepard plot shows scatter around the regression between the inter-point
#distances in the final configuration (i.e., the distances between each pair of communities)
#against their original dissimilarities
png("nmds_Species_Oct2024run.png", res = 300, width = 3600, height = 2700) #creates a file in your working folder to house the graph when it populates below
nmds.plot <- ordiplot(nmds2d,display="sites") #populates the plot into a variable
## Adding ellipses to group years
ordihull(nmds.plot,groups=met$SPP,draw="lines",col=c("tomato4","lightblue","springgreen3","black","blue3")) # adds ellipses around the point groups (OPTIONAL!)
##adjust colors to match each year, pch=20 makes it bullet points
#unique(met$SPP)
points(nmds.plot,"sites", pch=20, col= "tomato4", select = met$SPP == "BAST")     # set the metadata to the variable you need it to be
points(nmds.plot,"sites", pch=20, col= "lightblue", select = met$SPP == "ECAR")  # and set a different color to each level of the variable
points(nmds.plot,"sites", pch=20, col= "springgreen4", select = met$SPP == "EFLE")
points(nmds.plot,"sites", pch=20, col= "black", select = met$SPP == "MMUR")
points(nmds.plot,"sites", pch=20, col= "blue", select = met$SPP == "SW")
##Add Stress Value
text(0.56,0.35,"2D Stress: 0.207\nANOSIM R = 0.8608", cex=1.2) # make sure you add the stress value in an empty portion of the graph
##Adding legend
legend("topleft",legend= c("BAST","ECAR", "EFLE","MMUR","SW"),   # customize the legend to match the colors and variables
       title = "Species",
       col=c("tomato4","lightblue","springgreen4","black","blue3"),
       pch=19, cex=1.2)
##Adding title
title(main="nMDS of Relative Abundances by Octocoral Species and SW") # adds a title to the graph
graphics.off() #saves the graph to the file

# using qiime2R





head(taxonomy1) 

SVs<-read_qza("table-NoEMC.qza")$data
taxonomy<-read_qza("taxonomy.qza")$data %>% parse_taxonomy()
taxasums<-summarize_taxa(SVs, taxonomy)$Genus

met= readr::read_csv(file="metadata.csv")

view(taxasums)
taxa_barplot(taxasums, met, "SPP")
head(taxasums)

