## Root Growth Dynamics, Dominant Rhizosphere Bacteria, and Correlation Between Dominant Bacterial Genera and Root Traits Through Brassica Napus Development Taye et al. 2020

setwd("D:/.....")

##

library(ggplot2)
library(phyloseq)
library(ape)
library(biomformat)

#import tax, asv and metadata tables to create phyloseq object. Here the processed final asv and tax table is used. the only update is more sample vaiables (root traits) are added in the metadata

otu_table=read.csv("can.asv.rhizo_gg.csv", sep = ",",row.names = 1)
otu_table=as.matrix(otu_table)

taxonomy=read.csv("can.tax.rhizo_gg.csv",sep="," ,row.names=1)
taxonomy=as.matrix(taxonomy)

metadata=read.csv("can.sam.rhizo_gg.csv",sep="," ,row.names=1)

### import them as phyloseq objects

OTU=otu_table(otu_table,taxa_are_rows=TRUE)
TAX=tax_table(taxonomy)
META=sample_data(metadata)

## cheke that your otu names are consistent across objects

taxa_names(TAX)
taxa_names(OTU)

## make sure files have the same sample names

sample_names(OTU)
sample_names(META)


## Everything looks good --- now proceed to merging them into one phyloseq object

physeq.updated.line.names.new.variable=phyloseq(OTU,TAX,META)

sample_variables(physeq.updated.line.names.new.variable)## check if sample varibales are correctly presented in the                                                                    phyyloseq object

## save the updated phyloseq object

save(physeq.updated.line.names.new.variable,file="2016.physeq.updated.line.names.new.variable.RData")

### Changes in rhizosphere bacterial community associated to three different developmental stages of B. napus

setwd("D:.....")

library(microbiome)
library(ggpubr)
library(knitr)
library(dplyr)

### Dominant taxonomic groups per sample: at different canola lines, growth stages

# subset vegetative stage

veg = subset_samples(physeq.updated.line.names.new.variable, Growth.Stage=="Vegetative")
dp.1 = dominant(veg, level = "Phylum")
dp.2 = dominant(veg, level = "Class")
dp.3 = dominant(veg, level = "Order")
dp.4 = dominant(veg, level = "Family")
dp.5 = dominant(veg, level = "Genus")

## write dominant taxa at vegetative

write.csv(dp.1, 'veg.dominant.phylum.csv')
write.csv(dp.2, 'veg.dominant.class.csv')
write.csv(dp.3, 'veg.dominant.order.csv')
write.csv(dp.4, 'veg.dominant.family.csv')
write.csv(dp.5, 'veg.dominant.genus.csv')

#subset flowering

Flw = subset_samples(physeq.updated.line.names.new.variable, Growth.Stage=="Flowering")

dp.1 = dominant(Flw, level = "Phylum")
dp.2 = dominant(Flw, level = "Class")
dp.3 = dominant(Flw, level = "Order")
dp.4 = dominant(Flw, level = "Family")
dp.5 = dominant(Flw, level = "Genus")

## write dominant taxa at vegetative

write.csv(dp.1, 'Flw.dominant.phylum.csv')
write.csv(dp.2, 'Flw.dominant.class.csv')
write.csv(dp.3, 'Flw.dominant.order.csv')
write.csv(dp.4, 'Flw.dominant.family.csv')
write.csv(dp.5, 'Flw.dominant.genus.csv')

# subset maturity

mtu = subset_samples(physeq.updated.line.names.new.variable, Growth.Stage=="Maturity")

dp.1 = dominant(mtu, level = "Phylum")
dp.2 = dominant(mtu, level = "Class")
dp.3 = dominant(mtu, level = "Order")
dp.4 = dominant(mtu, level = "Family")
dp.5 = dominant(mtu, level = "Genus")

## write dominant taxa at vegetative
write.csv(dp.1, 'mtu.dominant.phylum.csv')
write.csv(dp.2, 'mtu.dominant.class.csv')
write.csv(dp.3, 'mtu.dominant.order.csv')
write.csv(dp.4, 'mtu.dominant.family.csv')
write.csv(dp.5, 'mtu.dominant.genus.csv')

## top ten taxa by mean abundance across growth stage

setwd("D:/out_puts/phyloseq/with_green_genes_taxa/Ampvis_visualization/Mean_abundance_growth_stage_ampvis")
setwd("D:/out_puts/phyloseq/with_green_genes_taxa/Ampvis_visualization")

library(ampvis2)
library(radiant) ## this is for rownames_to_column function
library(tibble)## for as.table function
library(phyloseq)

### make ampvis data object from phyloseq

ps_data = physeq.updated.line.names.new.variable ## my already created phyloseq object

## convert phyloseq object into ampvis object 
# Load data into ampvis

otu.table <- otu_table(ps_data)
tax.table <- tax_table(ps_data)
otu.table <- as.data.frame(otu.table)
tax.table <- as.data.frame(tax.table)
meta.table <- sample_data(ps_data)
meta.table <- as.tibble(meta.table)
write.csv(otu.table, 'otu.table.amp.csv')## and add OTU column name and reload
write.csv(tax.table, 'tax.table.amp.csv')

## load the OTUID name otu and tax files

otu.table.amp <- read.csv("otu.table.amp.csv")
tax.table.amp <- read.csv("tax.table.amp.csv")

#merge the two as datafrmae

otu.table.tax <- merge.data.frame(otu.table.amp, tax.table.amp, by= "OTU")
head(rownames(otu.table.tax))
head(rownames(otu.table.tax))

class(meta.table)
head(rownames(otu.table))
head(rownames(tax.table))
head(colnames(otu.table))
head(colnames(tax.table))

# Creat the Ampvis object

ps1.av2 <- amp_load(otutable = otu.table.tax, metadata = meta.table)

##cummulative rank abundance of genera at three growth stages

amp_rankabundance(ps1.av2, group_by = "Growth_Stage")

# Top ten most abudnant Phyla based on relative mean read abundance
## Phyla - across canola growth stages
# reordr the levels of growthsatge so that it goes from veg to mat in the plots

ps1.av2$metadata$Growth_Stage <- factor(ps1.av2$metadata$Growth_Stage, levels = c("Vegetative", "Flowering", "Maturity"))

# plot top 10 phylum across growth stages

amp_heatmap(ps1.av2, group_by = "Growth_Stage", plot_values = TRUE) +
  theme(axis.text.x = element_text(angle = 90, size=12, vjust = 1),
        axis.text.y = element_text(size=12),axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), 
        legend.position="right", legend.text = element_text(size = 12), legend.title = element_text(size = 12))## across                                                                                  canola growth stages## top 10 phyal

## Top 10 class

# plot top 10 class across growth stages
amp_heatmap(ps1.av2, group_by = "Growth_Stage", tax_aggregate = "Class", tax_add = "Phylum",plot_values = TRUE) +
  theme(axis.text.x = element_text(angle = 90, size=12, vjust = 1),axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size=12),
        legend.position="right", legend.text = element_text(size = 12), legend.title = element_text(size = 12))## across                                                                                                      canola growth stages

## Top 10 Family

# plot top 10 Family across growth stages

amp_heatmap(ps1.av2, group_by = "Growth_Stage", tax_aggregate = "Family", tax_add = "Phylum",plot_values = TRUE) +
  theme(axis.text.x = element_text(angle = 90, size=12, vjust = 1),axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size=12),
        legend.position="right", legend.text = element_text(size = 12), legend.title = element_text(size = 12))## across                                                                                                      canola growth stages

## Top 10 Genus

# plot top 10 genera across growth stages

amp_heatmap(ps1.av2, group_by = "Growth_Stage", tax_aggregate = "Genus", tax_add = "Phylum",plot_values = TRUE) +
  theme(axis.text.x = element_text(angle = 90, size=12, vjust = 1),
        axis.text.y = element_text(size=12),axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        legend.position="right", legend.text = element_text(size = 12), legend.title = element_text(size = 12))## across                                                                                                       canola growth stages

############## Frequency of each taxa at different growth stages #############

setwd("D:/out_puts/phyloseq/with_green_genes_taxa/Ampvis_visualization/Frequency_of_taxa_by_growth_stage")

# subset each growth stage

veg
Flw
mtu

# remove zero sum taxa

veg= prune_taxa(taxa_sums(veg) > 0, veg)
Flw= prune_taxa(taxa_sums(Flw) > 0, Flw)
mtu= prune_taxa(taxa_sums(mtu) > 0, mtu)

## frequecny  of phyla in each growth stage

num.taxa.py.vg = table(tax_table(veg)[,"Phylum"])
num.taxa.py.fl = table(tax_table(Flw)[,"Phylum"])
num.taxa.py.mt = table(tax_table(mtu)[,"Phylum"])

# create a dataframe

num.taxa.py.vg = as.data.frame(num.taxa.py.vg)
head(num.taxa.py.vg)

num.taxa.py.fl = as.data.frame(num.taxa.py.fl)
head(num.taxa.py.fl)

num.taxa.py.mt = as.data.frame(num.taxa.py.mt)
head(num.taxa.py.mt)

## save as csv

write.csv(num.taxa.py.vg, 'Phylum_frequency_vegeative.csv')
write.csv(num.taxa.py.fl, 'Phylum_frequency_flowering.csv')
write.csv(num.taxa.py.mt, 'Phylum_frequency_maturity.csv')

## Frequency of Class at each growth stage

num.taxa.cls.vg = table(tax_table(veg)[,"Class"])
num.taxa.cls.fl = table(tax_table(Flw)[,"Class"])
num.taxa.cls.mt = table(tax_table(mtu)[,"Class"])

# create a dataframe

num.taxa.cls.vg = as.data.frame(num.taxa.cls.vg)
head(num.taxa.cls.vg)

num.taxa.cls.fl = as.data.frame(num.taxa.cls.fl)
head(num.taxa.cls.fl)

num.taxa.cls.mt = as.data.frame(num.taxa.cls.mt)
head(num.taxa.cls.mt)

## save as csv

write.csv(num.taxa.cls.vg, 'Class_frequency_vegeative.csv')
write.csv(num.taxa.cls.fl, 'Class_frequency_flowering.csv')
write.csv(num.taxa.cls.mt, 'Class_frequency_maturity.csv')

## Frequence of order at eahc growth stage

num.taxa.ord.vg = table(tax_table(veg)[,"Order"])
num.taxa.ord.fl = table(tax_table(Flw)[,"Order"])
num.taxa.ord.mt = table(tax_table(mtu)[,"Order"])

# create a dataframe

num.taxa.ord.vg = as.data.frame(num.taxa.ord.vg)
head(num.taxa.ord.vg)

num.taxa.ord.fl = as.data.frame(num.taxa.ord.fl)
head(num.taxa.ord.fl)

num.taxa.ord.mt = as.data.frame(num.taxa.ord.mt)
head(num.taxa.ord.mt)

## save as csv

write.csv(num.taxa.ord.vg, 'Order_frequency_vegeative.csv')
write.csv(num.taxa.ord.fl, 'Order_frequency_flowering.csv')
write.csv(num.taxa.ord.mt, 'Order_frequency_maturity.csv')

## Frequency of Family at eahc growth stage

num.taxa.fam.vg = table(tax_table(veg)[,"Family"])
num.taxa.fam.fl = table(tax_table(Flw)[,"Family"])
num.taxa.fam.mt = table(tax_table(mtu)[,"Family"])

# create a dataframe

num.taxa.fam.vg = as.data.frame(num.taxa.fam.vg)
head(num.taxa.fam.vg)

num.taxa.fam.fl = as.data.frame(num.taxa.fam.fl)
head(num.taxa.fam.fl)

num.taxa.fam.mt = as.data.frame(num.taxa.fam.mt)
head(num.taxa.fam.mt)

## save as csv

write.csv(num.taxa.fam.vg, 'Family_frequency_vegeative.csv')
write.csv(num.taxa.fam.fl, 'Family_frequency_flowering.csv')
write.csv(num.taxa.fam.mt, 'Family_frequency_maturity.csv')

## Frequency of genus at eahc growth stage

num.taxa.gns.vg = table(tax_table(veg)[,"Genus"])
num.taxa.gns.fl = table(tax_table(Flw)[,"Genus"])
num.taxa.gns.mt = table(tax_table(mtu)[,"Genus"])

# create a dataframe

num.taxa.gns.vg = as.data.frame(num.taxa.gns.vg)
head(num.taxa.gns.vg)

num.taxa.gns.fl = as.data.frame(num.taxa.gns.fl)
head(num.taxa.gns.fl)

num.taxa.gns.mt = as.data.frame(num.taxa.gns.mt)
head(num.taxa.gns.mt)

## save as csv

write.csv(num.taxa.gns.vg, 'Genus_frequency_vegeative.csv')
write.csv(num.taxa.gns.fl, 'Genus_frequency_flowering.csv')
write.csv(num.taxa.gns.mt, 'Genus_frequency_maturity.csv')

##### ............................#### ...............##### ......................... #####

## Correlation between canola rhizosphere dominat genera and previously identified canola root-growth promoting 
#bacteria and canola root traits under field condition at different growth stages

setwd("D:/out_puts/phyloseq/with_green_genes_taxa/Ampvis_visualization/Root_traits_canola_PGPB")

## top ten genera in all growth stages based on percent of frequency of occurence

# Pseudomonas 
# Flavobacterium 
# Lysobacter 
# Rhodoplanes 
# Pedobacter 
# Vibrio 
# Stenotrophomonas 
# Arthrobacter 


## top ten genera in all growth stages based on percent mean relative abundance 

# Arthrobacter
# Stenotrophomonas
# Gluconacetobacter
# Skermanella
# Pseudomonas
# Bradyrhizobium 


## Total dominat genera (based on both frequence and mean relative abundance): "Use this Dominat genera for correlation"

# Pseudomonas 
# Flavobacterium 
# Lysobacter 
# Rhodoplanes 
# Pedobacter 
# Vibrio 
# Stenotrophomonas 
# Arthrobacter
# Gluconacetobacter
# Skermanella
# Bradyrhizobium 

## Canola PGPR previously indetified in canola (identified by Jorge et al 2019 and Bertrand et al., 2000): 
## " use this canola PGPR for coreelation" 

# Paenibacillus 
# Pantoea 
# Pseudomonas 
# Rhizobium 
# Stenotrophomonas 
# agrococcus 
# Bacillus 
# Leifsonia 
# Microbacterium 
# Rhodococcus 
# Xanthomonas
# Variovorax 
# Agrobacterium 
# Phyllobacterium

####load libraries 

library(ggplot2)
library(phyloseq)
library(ape)
library(biomformat)
library(microbiome)
library(ggpubr)
library(knitr)
library(dplyr)
library(ampvis2)
library(radiant) ## this is for rownames_to_column function
library(tibble)## for as.table function
library(microbiomeSeq)


# load the phyloseq object containing the  2016 rhizosphere bacterial dataset

load("D:/out_puts/phyloseq/with_green_genes_taxa/2016.physeq.updated.line.names.new.variable.RData")

#### subset dominat genera based on both frequence and mean relative abundance


physeq.top.10.genera.root.trait =  subset_taxa(physeq.updated.line.names.new.variable, Genus=="Pseudomonas"| Genus=="Flavobacterium"|
                                                 Genus=="Lysobacter"|Genus=="Rhodoplanes"|Genus=="Pedobacter"|Genus=="Vibrio"|Genus=="Stenotrophomonas"
                                               |Genus=="Arthrobacter"|Genus=="Gluconacetobacter"|Genus=="Skermanella"|Genus=="Bradyrhizobium")

## remove zero sum taxa

physeq.top.10.genera.root.trait= prune_taxa(taxa_sums(physeq.top.10.genera.root.trait) > 0, physeq.top.10.genera.root.trait)

## Using the taxa.env.correlation function of microbiome seq, lets do correlation between genera and root traits

sample_variables(physeq.top.10.genera.root.trait) ## get the sample variable names

physeq <- taxa_level(physeq.top.10.genera.root.trait, "Genus") ### tax_glom at geneus level with genus name indicated


# do the correlation at each growth stage (veg, flw and Mat)

root.top.10.taxa.cor <- taxa.env.correlation(physeq, grouping_column = "Growth.Stage", method = "pearson", 
                                             pvalue.threshold = 0.05, padjust.method = "BH", adjustment = 4, num.taxa = 173, 
                                             select.variables = c("PercentExtraFRL", "PercentFRL" , "Root.Biomass", "Total_RL","CoarseRL" ,"TFRL2mm", "ExtraFRL", "FRL"))


# reorder growth stage levels

root.top.10.taxa.cor$Type <- factor(root.top.10.taxa.cor$Type, levels = c("Vegetative", "Flowering", "Maturity"))

p <- plot_taxa_env(root.top.10.taxa.cor)
p = p + xlab("Growth stages")
p = p + theme(axis.title.x = element_text(size = 12, face = "bold"), axis.text.x=element_text(size = 12))
p = p + theme(axis.text.y=element_text(size = 12))
p = p + theme(axis.text.y=element_text(angle = 45))

print(p) #  p-values for multiple comparisions using Benjamin and Hochberg: adjustment = 4 - adjust Taxa (row on the                  correlation plot)

#  save csv file containing the correlation and p values

write.csv(root.top.10.taxa.cor, 'dominat.genera.root.trait.corr.with.pvalues.csv')

### Subset Canola PGPR previously indetified in canola (identified by Jorge et al 2019 and Bertrand et al., 2000)
## of the 14 pgpr reported 12 of thee genera were present in our dataset 

physeq.pgpr.genera.root.trait =  subset_taxa(physeq.updated.line.names.new.variable, Genus=="Paenibacillus"| Genus=="Pantoea"|
                                               Genus=="Pseudomonas"|Genus=="Rhizobium"|Genus=="Stenotrophomonas"|Genus=="agrococcus"|Genus=="Bacillus"
                                             |Genus=="Leifsonia"|Genus=="Microbacterium"|Genus=="Rhodococcus"|Genus=="Xanthomonas"|Genus=="Variovorax"
                                             |Genus=="Agrobacterium"|Genus=="Phyllobacterium")

## remove zero sum taxa

physeq.pgpr.genera.root.trait= prune_taxa(taxa_sums(physeq.pgpr.genera.root.trait) > 0, physeq.pgpr.genera.root.trait)

## Using the taxa.env.correlation function of microbiome seq, lets do correlation between genera and root traits

sample_variables(physeq.pgpr.genera.root.trait) ## get the sample variable names

physeq1 <- taxa_level(physeq.pgpr.genera.root.trait, "Genus") ### tax_glom at geneus level with genus name indicated


# do the correlation at each growth stage (veg, flw and Mat)

root.pgpr.taxa.cor <- taxa.env.correlation(physeq1, grouping_column = "Growth.Stage", method = "pearson", 
                                           pvalue.threshold = 0.05, padjust.method = "BH", adjustment = 4, num.taxa = 173, 
                                           select.variables = c("PercentExtraFRL", "PercentFRL" , "Root.Biomass", "Total_RL","CoarseRL" ,"TFRL2mm", "ExtraFRL", "FRL"))


# reorder growth stage levels

root.pgpr.taxa.cor$Type <- factor(root.pgpr.taxa.cor$Type, levels = c("Vegetative", "Flowering", "Maturity"))

p1 <- plot_taxa_env(root.pgpr.taxa.cor)
p1 = p1 + xlab("Growth stages")
p1 = p1 + theme(axis.title.x = element_text(size = 12, face = "bold"), axis.text.x=element_text(size = 12))
p1 = p1 + theme(axis.text.y=element_text(size = 12))
p1 = p1 + theme(axis.text.y=element_text(angle = 45))

print(p1) #  p-values for multiple comparisions using Benjamin and Hochberg: adjustment = 4 - adjust Taxa (row on the correlation plot)

#  save csv file containing the correlation and p values

write.csv(root.pgpr.taxa.cor, 'pgprt.genera.root.trait.corr.with.pvalues.csv')

### Correlation of dominat and canola gprb with environmental and soil characteristics


physeq.top.10.genera.root.trait =  subset_taxa(physeq.updated.line.names.new.variable, Genus=="Pseudomonas"| Genus=="Flavobacterium"|
                                                 Genus=="Lysobacter"|Genus=="Rhodoplanes"|Genus=="Pedobacter"|Genus=="Vibrio"|Genus=="Stenotrophomonas"
                                               |Genus=="Arthrobacter"|Genus=="Gluconacetobacter"|Genus=="Skermanella"|Genus=="Bradyrhizobium")

## remove zero sum taxa

physeq.top.10.genera.root.trait= prune_taxa(taxa_sums(physeq.top.10.genera.root.trait) > 0, physeq.top.10.genera.root.trait)

physeq <- taxa_level(physeq.top.10.genera.root.trait, "Genus") ### tax_glom at geneus level with genus name indicated

# do the correlation at each growth stage (veg, flw and Mat)

envi.top.10.taxa.cor <- taxa.env.correlation(physeq, grouping_column = "Growth.Stage", method = "pearson", 
                                             pvalue.threshold = 0.05, padjust.method = "BH", adjustment = 4, num.taxa = 173, 
                                             select.variables = c("T.WP", "PPT.WP", "T.SD","PPT.SD" ,"Soil.moisture", "Soil.water.content"))


# reorder growth stage levels
envi.top.10.taxa.cor$Type <- factor(envi.top.10.taxa.cor$Type, levels = c("Vegetative", "Flowering", "Maturity"))

p2 <- plot_taxa_env(envi.top.10.taxa.cor)
p2 = p2 + xlab("Growth stages")
p2 = p2 + theme(axis.title.x = element_text(size = 12, face = "bold"), axis.text.x=element_text(size = 12))
p2 = p2 + theme(axis.text.y=element_text(size = 12))
p2 = p2 + theme(axis.text.y=element_text(angle = 45))

print(p2) #  p-values for multiple comparisions using Benjamin and Hochberg: adjustment = 4 - adjust Taxa (row on the                  correlation plot)

#  save csv file containing the correlation and p values

write.csv(envi.top.10.taxa.cor, 'ominat.genera.envi.soil.variable.corr.with.pvalues.csv')

##

physeq.pgpr.genera.root.trait =  subset_taxa(physeq.updated.line.names.new.variable, Genus=="Paenibacillus"| Genus=="Pantoea"|
                                               Genus=="Pseudomonas"|Genus=="Rhizobium"|Genus=="Stenotrophomonas"|Genus=="agrococcus"|Genus=="Bacillus"
                                             |Genus=="Leifsonia"|Genus=="Microbacterium"|Genus=="Rhodococcus"|Genus=="Xanthomonas"|Genus=="Variovorax"
                                             |Genus=="Agrobacterium"|Genus=="Phyllobacterium")

## remove zero sum taxa

physeq.pgpr.genera.root.trait= prune_taxa(taxa_sums(physeq.pgpr.genera.root.trait) > 0, physeq.pgpr.genera.root.trait)

## Using the taxa.env.correlation function of microbiome seq, lets do correlation between genera and root traits

sample_variables(physeq.pgpr.genera.root.trait) ## get the sample variable names

physeq1 <- taxa_level(physeq.pgpr.genera.root.trait, "Genus") ### tax_glom at geneus level with genus name indicated


# do the correlation at each growth stage (veg, flw and Mat)

envi.pgpr.taxa.cor <- taxa.env.correlation(physeq1, grouping_column = "Growth.Stage", method = "pearson", 
                                           pvalue.threshold = 0.05, padjust.method = "BH", adjustment = 4, num.taxa = 173, 
                                           select.variables = c("T.WP", "PPT.WP", "T.SD","PPT.SD" ,"Soil.moisture", "Soil.water.content"))


# reorder growth stage levels

envi.pgpr.taxa.cor$Type <- factor(envi.pgpr.taxa.cor$Type, levels = c("Vegetative", "Flowering", "Maturity"))

p3 <- plot_taxa_env(envi.pgpr.taxa.cor)
p3 = p3 + xlab("Growth stages")
p3 = p3 + theme(axis.title.x = element_text(size = 12, face = "bold"), axis.text.x=element_text(size = 12))
p3 = p3 + theme(axis.text.y=element_text(size = 12))
p3 = p3 + theme(axis.text.y=element_text(angle = 45))

print(p3) #  p-values for multiple comparisions using Benjamin and Hochberg: adjustment = 4 - adjust Taxa (row on the                 correlation plot)

#  save csv file containing the correlation and p values

write.csv(envi.pgpr.taxa.cor, 'pgprt.genera.envi.soil.variables.corr.with.pvalues.csv')


#### ....... ####


### Ratio among the three dominat phyla: how the ratio changes across growth stages + correlation between acidobacteria:proteobactera to root traits

## Agglomerate taxa to phylum level

Physeq.phylum = tax_glom(physeq.updated.line.names.new.variable, taxrank = "Phylum")

##subset the two most abundant phyla: Acidobacteria and proteobacteria

Physeq.phylum.aciso.proteo = subset_taxa(Physeq.phylum, Phylum=="Acidobacteria"|Phylum=="Proteobacteria")

# get abudnace of Acidobacteria at the three growth stages

## modify the source microbiomer r function (brratio) that calculate the ratio Bacteriodetes/Firmicutes

library(microbiome)

bfratio <- function (x) {
  a <- transform(aggregate_taxa(x, level = "Phylum"), "compositional")
  b <- abundances(a)["Proteobacteria",]
  f <- abundances(a)["Acidobacteria",]  
  #data.frame(Bacteroidetes = b, Firmicutes = f, Ratio = b/f)
  b/f
}

## subset vegetative, flowering and maturity stage samples

Physeq.phylum.aciso.proteo.veg = subset_samples(Physeq.phylum.aciso.proteo, Growth.Stage=="Vegetative")
Physeq.phylum.aciso.proteo.flw = subset_samples(Physeq.phylum.aciso.proteo, Growth.Stage=="Flowering")
Physeq.phylum.aciso.proteo.mat = subset_samples(Physeq.phylum.aciso.proteo, Growth.Stage=="Maturity")

## calculate Proteobacteria/Acidobacteria ratio

mean(bfratio(Physeq.phylum.aciso.proteo.veg)) #  5.094917

pro.acid.ratio= bfratio(Physeq.phylum.aciso.proteo.flw) # there is one sample (ORIG1CS1143NW06000c) with inf value (no                                                               acidobacteria). Remove it to calculate mean

pro.acid.ratio = pro.acid.ratio[!is.infinite(pro.acid.ratio)]
mean(pro.acid.ratio)# flowering sateg P/A ratio =  10.26853

mean(bfratio(Physeq.phylum.aciso.proteo.mat)) # 7.939675

## save proteobacteria:Acidobacteria as a datarame and correlate with fine root traits.

veg.Pro.acido.ratio = bfratio(Physeq.phylum.aciso.proteo.veg)

write.csv(veg.Pro.acido.ratio, 'veg.Pro.acido.ratio.csv')## add column nammes and re-read to merge with metadata
veg.Pro.acido.ratio= read.csv("veg.Pro.acido.ratio.csv")
veg.Pro.acido.ratio= as.data.frame(veg.Pro.acido.ratio)
metav = sample_data(Physeq.phylum.aciso.proteo.veg)
metav = as.data.frame(metav)
metav$P.A.ratio = veg.Pro.acido.ratio$P.A.Ratio ## add the P:A ratio to the meta data table containing the root traits.

#
flw.pro.acid.ratio= bfratio(Physeq.phylum.aciso.proteo.flw)
write.csv(flw.pro.acid.ratio, 'flw.Pro.acido.ratio.csv')## add column nammes and re-read to merge with metadata
flw.pro.acid.ratio= read.csv("flw.Pro.acido.ratio.csv")
flw.pro.acid.ratio= as.data.frame(flw.pro.acid.ratio)
metaf = sample_data(Physeq.phylum.aciso.proteo.flw)
metaf = as.data.frame(metaf)
metaf$P.A.ratio = flw.pro.acid.ratio$P.A.Ratio ## add the P:A ratio to the meta data table containing the root traits.

##
mat.pro.acid.ratio= bfratio(Physeq.phylum.aciso.proteo.mat)
write.csv(mat.pro.acid.ratio, 'mat.pro.acid.ratio.csv')## add column nammes and re-read to merge with metadata
mat.pro.acid.ratio= read.csv("mat.pro.acid.ratio.csv")
mat.pro.acid.ratio= as.data.frame(mat.pro.acid.ratio)
metam = sample_data(Physeq.phylum.aciso.proteo.mat)
metam = as.data.frame(metam)
metam$P.A.ratio = mat.pro.acid.ratio$P.A.Ratio ## add the P:A ratio to the meta data table containing the root traits.

#### plot correlation between root traits and P:A ratio

library(ggpubr)

p= ggscatter(metav, x = "Total_RL", y = "P.A.ratio", 
             add = "reg.line", conf.int = TRUE,
             cor.coef = TRUE, cor.method = "pearson",
             xlab="Total root length", ylab = "Proteobacteria:Acidobacteria")## scater plot with trend line and correlation                                                                                 coefficient with p value 
p= p +theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14))
p=p+theme(axis.text.x = element_text(size=12), axis.text.y=element_text(size = 12))
p


##

p1= ggscatter(metav, x = "CoarseRL", y = "P.A.ratio", 
              add = "reg.line", conf.int = TRUE,
              cor.coef = TRUE, cor.method = "pearson",
              xlab="Coarse root length", ylab = "")## scater plot with trend line and correlation coefficient with p value 
p1= p1 +theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14))
p1=p1+theme(axis.text.x = element_text(size=12), axis.text.y=element_text(size = 12))
p1

##

p2= ggscatter(metav, x = "TFRL2mm", y = "P.A.ratio", 
              add = "reg.line", conf.int = TRUE,
              cor.coef = TRUE, cor.method = "pearson",
              xlab="Total Fine root length", ylab = "Proteobacteria:Acidobacteria")## scater plot with trend line and                                                                                          correlation coefficient with p value 
p2= p2 +theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14))
p2=p2+theme(axis.text.x = element_text(size=12), axis.text.y=element_text(size = 12))
p2

##

p3= ggscatter(metav, x = "Root.Biomass", y = "P.A.ratio", 
              add = "reg.line", conf.int = TRUE,
              cor.coef = TRUE, cor.method = "pearson",
              xlab="Root biomass", ylab = "")## scater plot with trend line and correlation coefficient with p value 
p3= p3 +theme(axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14))
p3=p3+theme(axis.text.x = element_text(size=12), axis.text.y=element_text(size = 12))
p3

###
ggarrange(p, p1, p2, p3, labels = c("A", "B","C", "D"))

#### Actinobacteria: Acidobacteria ratio

aaratio <- function (x) {
  a <- transform(aggregate_taxa(x, level = "Phylum"), "compositional")
  b <- abundances(a)["Actinobacteria",]
  f <- abundances(a)["Acidobacteria",]  
  #data.frame(Bacteroidetes = b, Firmicutes = f, Ratio = b/f)
  b/f
}

##

Physeq.phylum.acti.acido = subset_taxa(Physeq.phylum, Phylum=="Acidobacteria"|Phylum=="Actinobacteria")

## subset vegetative, flowering and maturity stage samples
Physeq.phylum.acti.acido.veg = subset_samples(Physeq.phylum.acti.acido, Growth.Stage=="Vegetative")
Physeq.phylum.acti.acido.flw = subset_samples(Physeq.phylum.acti.acido, Growth.Stage=="Flowering")
Physeq.phylum.acti.acido.mat = subset_samples(Physeq.phylum.acti.acido, Growth.Stage=="Maturity")
##

mean(aaratio(Physeq.phylum.acti.acido.veg))#2.027333
acti.acid.ratio = aaratio(Physeq.phylum.acti.acido.flw)# there is one sample (ORIG1CS1143NW06000c) with inf value (no acidobacteria). Remove it to calculate mean
acti.acid.ratio = acti.acid.ratio[!is.infinite(acti.acid.ratio)]
mean(acti.acid.ratio)#3.530675 #flowering

mean(aaratio(Physeq.phylum.acti.acido.mat))# 4.224899

############## genus level P:A RATIO ################################

Physeq.class = tax_glom(physeq.updated.line.names.new.variable, taxrank = "Genus")

##subset the two most abundant phyla: Acidobacteria and proteobacteria

Physeq.steno_arthro = subset_taxa(Physeq.class, Genus=="Stenotrophomonas"|Genus=="Arthrobacter")
Physeq.skerm_arthro = subset_taxa(Physeq.class, Genus=="Skermanella"|Genus=="Arthrobacter")
Physeq.Pseudo_arthro = subset_taxa(Physeq.class, Genus=="Bradyrhizobium"|Genus=="Arthrobacter")


#get_taxa_unique(Physeq.class.aciso.proteo, taxonomic.rank=rank_names(Physeq.class.aciso.proteo)[3], errorIfNULL=TRUE)

# get abudnace of Acidobacteria at the three growth stages

## modify the source microbiomer r function (brratio) that calculate the ratio Bacteriodetes/Firmicutes
library(microbiome)
bfratio <- function (x) {
  a <- transform(aggregate_taxa(x, level = "Genus"), "compositional")
  b <- abundances(a)["Bradyrhizobium",]
  f <- abundances(a)["Arthrobacter",]  
  #data.frame(Bacteroidetes = b, Firmicutes = f, Ratio = b/f)
  b/f
}
## subset vegetative, flowering and maturity stage samples
Physeq.steno_arthro.veg = subset_samples(Physeq.steno_arthro, Growth.Stage=="Vegetative")
Physeq.steno_arthro.flw = subset_samples(Physeq.steno_arthro, Growth.Stage=="Flowering")
Physeq.steno_arthro.mat = subset_samples(Physeq.steno_arthro, Growth.Stage=="Maturity")

## calculate Proteobacteria/Acidobacteria ratio
alphapro.steno_arthro.ratio= bfratio(Physeq.steno_arthro.veg)
alphapro.steno_arthro.ratio = alphapro.steno_arthro.ratio[!is.nan(alphapro.steno_arthro.ratio)]
alphapro.steno_arthro.ratio = alphapro.steno_arthro.ratio[!is.infinite(alphapro.steno_arthro.ratio)]
mean(alphapro.steno_arthro.ratio)# vegetative sateg steno/Arthro ratio =  0.4964963

alphapro.steno_arthro.ratio= bfratio(Physeq.steno_arthro.flw) # 
alphapro.steno_arthro.ratio = alphapro.steno_arthro.ratio[!is.nan(alphapro.steno_arthro.ratio)]
alphapro.steno_arthro.ratio = alphapro.steno_arthro.ratio[!is.infinite(alphapro.steno_arthro.ratio)]
mean(alphapro.steno_arthro.ratio)# flowering sateg steno/Arthro ratio =  0.9993716

alphapro.steno_arthro.ratio= bfratio(Physeq.steno_arthro.mat) # 
alphapro.steno_arthro.ratio = alphapro.steno_arthro.ratio[!is.nan(alphapro.steno_arthro.ratio)]
alphapro.steno_arthro.ratio = alphapro.steno_arthro.ratio[!is.infinite(alphapro.steno_arthro.ratio)]
mean(alphapro.steno_arthro.ratio)# Maturity sateg steno/Arthro ratio =  0.5879472

##

## subset vegetative, flowering and maturity stage samples
Physeq.steno_arthro.veg = subset_samples(Physeq.skerm_arthro, Growth.Stage=="Vegetative")
Physeq.steno_arthro.flw = subset_samples(Physeq.skerm_arthro, Growth.Stage=="Flowering")
Physeq.steno_arthro.mat = subset_samples(Physeq.skerm_arthro, Growth.Stage=="Maturity")

## calculate Proteobacteria/Acidobacteria ratio
alphapro.steno_arthro.ratio= bfratio(Physeq.steno_arthro.veg)
alphapro.steno_arthro.ratio = alphapro.steno_arthro.ratio[!is.nan(alphapro.steno_arthro.ratio)]
alphapro.steno_arthro.ratio = alphapro.steno_arthro.ratio[!is.infinite(alphapro.steno_arthro.ratio)]
mean(alphapro.steno_arthro.ratio)# vegetative sateg Skerm/Arthro ratio =  1.152353

alphapro.steno_arthro.ratio= bfratio(Physeq.steno_arthro.flw) # 
alphapro.steno_arthro.ratio = alphapro.steno_arthro.ratio[!is.nan(alphapro.steno_arthro.ratio)]
alphapro.steno_arthro.ratio = alphapro.steno_arthro.ratio[!is.infinite(alphapro.steno_arthro.ratio)]
mean(alphapro.steno_arthro.ratio)# flowering sateg Skerm/Arthro ratio =  0.4422514

alphapro.steno_arthro.ratio= bfratio(Physeq.steno_arthro.mat) # 
alphapro.steno_arthro.ratio = alphapro.steno_arthro.ratio[!is.nan(alphapro.steno_arthro.ratio)]
alphapro.steno_arthro.ratio = alphapro.steno_arthro.ratio[!is.infinite(alphapro.steno_arthro.ratio)]
mean(alphapro.steno_arthro.ratio)# Maturity sateg steno/Arthro ratio =  0.4416642

##

## subset vegetative, flowering and maturity stage samples
Physeq.steno_arthro.veg = subset_samples(Physeq.Pseudo_arthro, Growth.Stage=="Vegetative")
Physeq.steno_arthro.flw = subset_samples(Physeq.Pseudo_arthro, Growth.Stage=="Flowering")
Physeq.steno_arthro.mat = subset_samples(Physeq.Pseudo_arthro, Growth.Stage=="Maturity")

## calculate Proteobacteria/Acidobacteria ratio
alphapro.steno_arthro.ratio= bfratio(Physeq.steno_arthro.veg)
alphapro.steno_arthro.ratio = alphapro.steno_arthro.ratio[!is.nan(alphapro.steno_arthro.ratio)]
alphapro.steno_arthro.ratio = alphapro.steno_arthro.ratio[!is.infinite(alphapro.steno_arthro.ratio)]
mean(alphapro.steno_arthro.ratio)# vegetative sateg Brady/Arthro ratio =  0.9194767

alphapro.steno_arthro.ratio= bfratio(Physeq.steno_arthro.flw) # 
alphapro.steno_arthro.ratio = alphapro.steno_arthro.ratio[!is.nan(alphapro.steno_arthro.ratio)]
alphapro.steno_arthro.ratio = alphapro.steno_arthro.ratio[!is.infinite(alphapro.steno_arthro.ratio)]
mean(alphapro.steno_arthro.ratio)# flowering sateg Brady/Arthro ratio =  0.3059721

alphapro.steno_arthro.ratio= bfratio(Physeq.steno_arthro.mat) # 
alphapro.steno_arthro.ratio = alphapro.steno_arthro.ratio[!is.nan(alphapro.steno_arthro.ratio)]
alphapro.steno_arthro.ratio = alphapro.steno_arthro.ratio[!is.infinite(alphapro.steno_arthro.ratio)]
mean(alphapro.steno_arthro.ratio)# Maturity sateg Brady/Arthro ratio =  0.2874067


#### ........................................ #### ...................... #### ..................... ####

## Root growth dynamics 

setwd("D:/out_puts/phyloseq/with_green_genes_taxa/Ampvis_visualization/Root_traits_canola_PGPB")

#load the meta data with root traits

head(can_sam_rhizo_gg)
data.root = as.data.frame(can.sam.rhizo_gg)

##
#hist(data.root$Total_RL)
#hist(sqrt(data.root$FRL)) ## square root transformation normalize the data
#hist(sqrt(data.root$TFRL2mm)) 
#hist(sqrt(data.root$CoarseRL)) 
#hist(sqrt(data.root$ExtraFRL)) 
#hist(sqrt(data.root$Total_RL)) 

#### visualize the data: by plotting
#with(data.root, cor(FRL, TFRL2mm))
#ggplot(data=data.root, aes(x=Growth.Stage, y=FRL))+ geom_point()

## fit generalized additive models for each root lenght trait to evaluate root dynaics throughout the ten sampling weeks
##fine root length:FRL
## GAM model using whole data and sampling point
# fine root length
#Growth= as.factor(data.root$Growth.Stage)

library(mgcv)

#library(gam)
#week=as.factor(data.root$Week)

gm1 <- gam(FRL ~ s(Week),data=data.root)
plot(gm1)
summary(gm1)

pd_0 <- data.frame(Week=seq(1,10,length=477))

pv_0 <- predict(gm1,newdata=pd_0,type="link",se=TRUE)
plot(data.root$Week,data.root$FRL,type="p",xlab="Sampling Week",ylab="Total Fine Root Length (cm)",cex.axis=1.5,cex.lab=1.5, lwd=1)
lines((pd_0[,1]),pv_0$fit,lwd=3,col="blue")
lines((pd_0[,1]),pv_0$fit+pv_0$se.fit*1.96,lwd=1,lty=3)
lines((pd_0[,1]),pv_0$fit-pv_0$se.fit *1.96,lwd=1,lty=3)

plt= list(data.root$Plot)
gmx <- gamm(FRL ~ s(Week),random=list(plt=~1),data=data.root)

##TFRL2mm : fine root length 0.5 to 2mm

gm2 <- gam(TFRL2mm ~ s(Week)+s(plot,bs="re"),data=data.root)


gm2 <- gam(TFRL2mm ~ s(Week),data=data.root)
plot(gm2)
summary(gm2)

pd_1 <- data.frame(Week=seq(1,10,length=477))

pv_1 <- predict(gm2,newdata=pd_1,type="link",se=TRUE)
plot(data.root$Week,data.root$TFRL2mm,type="p",xlab="Sampling Week",ylab="Fine Root Length (0.5-2 mm) (cm)",cex.axis=1.5,cex.lab=1.5, lwd=1)
lines((pd_1[,1]),pv_1$fit,lwd=3,col="blue")
lines((pd_1[,1]),pv_1$fit+pv_0$se.fit*1.96,lwd=1,lty=3)
lines((pd_1[,1]),pv_1$fit-pv_0$se.fit *1.96,lwd=1,lty=3)

### ExtraFRL : extra fine root length (0 -0.5mm)

gm3 <- gam(ExtraFRL ~ s(Week),data=data.root)
plot(gm3)
summary(gm3)

pd_2 <- data.frame(Week=seq(1,10,length=477))

pv_2 <- predict(gm3,newdata=pd_2,type="link",se=TRUE)
plot(data.root$Week,data.root$ExtraFRL,type="p",xlab="Sampling Week",ylab="Extra Fine Root Length (0-0.5 mm) (cm)",cex.axis=1.5,cex.lab=1.5, lwd=1)
lines((pd_2[,1]),pv_2$fit,lwd=3,col="blue")
lines((pd_2[,1]),pv_2$fit+pv_0$se.fit*1.96,lwd=1,lty=3)
lines((pd_2[,1]),pv_2$fit-pv_0$se.fit *1.96,lwd=1,lty=3)

##CoarseRL : coarse rot length (cm)

gm4 <- gam(CoarseRL ~ s(Week),data=data.root)
plot(gm4)
summary(gm4)

pd_3 <- data.frame(Week=seq(1,10,length=477))

pv_3 <- predict(gm4,newdata=pd_3,type="link",se=TRUE)
plot(data.root$Week,data.root$CoarseRL,type="p",xlab="Sampling Week",ylab="Coarse root length (cm)",cex.axis=1.5,cex.lab=1.5, lwd=1)
lines((pd_3[,1]),pv_3$fit,lwd=3,col="blue")
lines((pd_3[,1]),pv_3$fit+pv_0$se.fit*1.96,lwd=1,lty=3)
lines((pd_3[,1]),pv_3$fit-pv_0$se.fit *1.96,lwd=1,lty=3)

##Total_RL: 

gm5 <- gam(Total_RL ~ s(Week),data=data.root)
plot(gm5)
summary(gm5)

pd_4 <- data.frame(Week=seq(1,10,length=477))

pv_4 <- predict(gm5,newdata=pd_4,type="link",se=TRUE)
plot(data.root$Week,data.root$Total_RL,type="p",xlab="Sampling Week",ylab="Total root length (cm)",cex.axis=1.5,cex.lab=1.5, lwd=1)
lines((pd_4[,1]),pv_4$fit,lwd=3,col="blue")
lines((pd_4[,1]),pv_4$fit+pv_0$se.fit*1.96,lwd=1,lty=3)
lines((pd_4[,1]),pv_4$fit-pv_0$se.fit *1.96,lwd=1,lty=3)

##Root biomass

rootbiomass= data.root$`Root Biomass`
gm6 <- gam(rootbiomass ~ s(Week),data=data.root)
plot(gm6)
summary(gm6)

pd_5 <- data.frame(Week=seq(1,10,length=477))

pv_5 <- predict(gm6,newdata=pd_5,type="link",se=TRUE)
plot(data.root$Week,data.root$`Root Biomass`,type="p",xlab="Sampling Week",ylab="Root Biomass (g)",cex.axis=1.5,cex.lab=1.5, lwd=1)
lines((pd_5[,1]),pv_5$fit,lwd=3,col="blue")
lines((pd_5[,1]),pv_5$fit+pv_0$se.fit*1.96,lwd=1,lty=3)
lines((pd_5[,1]),pv_5$fit-pv_0$se.fit *1.96,lwd=1,lty=3)

#### plot together 

par(mfrow = c(2,2),xpd = NA,
    oma = c(4,4,2,0),
    mar = c(2,3,2,3))

plot(data.root$Week,data.root$FRL,type="p",xlab="",ylab ="Total Fine Root Length (cm)",cex.axis=1.5,cex.lab=1.5, lwd=1)
p+lines((pd_0[,1]),pv_0$fit,lwd=3,col="blue")
lines((pd_0[,1]),pv_0$fit+pv_0$se.fit*1.96,lwd=1,lty=3)
lines((pd_0[,1]),pv_0$fit-pv_0$se.fit *1.96,lwd=1,lty=3)
text(1,180,"(a)",cex=1.5)

plot(data.root$Week,data.root$TFRL2mm,type="p",xlab="",ylab="Fine Root Length (0.5-2 mm) (cm)",cex.axis=1.5,cex.lab=1.5, lwd=1)
lines((pd_1[,1]),pv_1$fit,lwd=3,col="blue")
lines((pd_1[,1]),pv_1$fit+pv_1$se.fit*1.96,lwd=1,lty=3)
lines((pd_1[,1]),pv_1$fit-pv_1$se.fit *1.96,lwd=1,lty=3)
text(1,1000,"(b)",cex=1.5)

plot(data.root$Week,data.root$ExtraFRL,type="p",xlab="Sampling Week",ylab="Extra Fine Root Length (0-0.5 mm) (cm)",cex.axis=1.5,cex.lab=1.5, lwd=1)
lines((pd_2[,1]),pv_2$fit,lwd=3,col="blue")
lines((pd_2[,1]),pv_2$fit+pv_2$se.fit*1.96,lwd=1,lty=3)
lines((pd_2[,1]),pv_2$fit-pv_2$se.fit *1.96,lwd=1,lty=3)
text(1,850,"(c)",cex=1.5)

plot(data.root$Week,data.root$CoarseRL,type="p",xlab="Sampling Week",ylab="Coarse root length (cm)",cex.axis=1.5,cex.lab=1.5, lwd=1)
lines((pd_3[,1]),pv_3$fit,lwd=3,col="blue")
lines((pd_3[,1]),pv_3$fit+pv_3$se.fit*1.96,lwd=1,lty=3)
lines((pd_3[,1]),pv_3$fit-pv_3$se.fit *1.96,lwd=1,lty=3)
text(1,50,"(d)",cex=1.5)

par(mfrow = c(1,2),xpd = NA,
    oma = c(4,4,2,0),
    mar = c(2,3,2,3))
plot(data.root$Week,data.root$Total_RL,type="p",xlab="Sampling Week",ylab="Total root length (cm)",cex.axis=1.5,cex.lab=1.5, lwd=1)
lines((pd_4[,1]),pv_4$fit,lwd=3,col="blue")
lines((pd_4[,1]),pv_4$fit+pv_4$se.fit*1.96,lwd=1,lty=3)
lines((pd_4[,1]),pv_4$fit-pv_4$se.fit *1.96,lwd=1,lty=3)
text(1,1000,"(a)",cex=1.5)

plot(data.root$Week,data.root$`Root Biomass`,type="p",xlab="Sampling Week",ylab="Root Biomass (g)",cex.axis=1.5,cex.lab=1.5, lwd=1)
lines((pd_5[,1]),pv_5$fit,lwd=3,col="blue")
lines((pd_5[,1]),pv_5$fit+pv_5$se.fit*1.96,lwd=1,lty=3)
lines((pd_5[,1]),pv_5$fit-pv_5$se.fit *1.96,lwd=1,lty=3)
text(1,4,"(b)",cex=1.5)

### Is there variability in root growth dynamics among genotypes. gam models considering canola genotypes: do on total root length and biomass

canola= as.factor(data.root$Canola.Lines)
plot= as.factor(data.root$Plot)

# fit gam model 

gm6 <- gam(Total_RL ~ s(Week, by = canola),data=data.root)
summary(gm6)
plot(gm6)


gm6 <- gam(FRL ~ s(Week,by = canola ),data=data.root)
plot(gm6)
summary(gm6)

par(mfrow = c(2,2),
    oma = c(5,3,2,0),
    mar = c(4,4,2,3))

pd_6 <- data.frame(Week=seq(1,10,length=30),canola=rep("NAM-0",30))
pv_6 <- predict(gm6,newdata=pd_6,type="link",se=TRUE)
plot(data.root$Week,data.root$Total_RL,type="p",xlab="",ylab="Total Root Length (cm)",cex.axis=1.5,cex.lab=1.5, lwd=1)
with(data.root[canola=="NAM-0",],points(Week,subset(Total_RL, canola=="NAM-0")))
lines((pd_6[,1]),pv_6$fit,lwd=3,col="blue")
lines((pd_6[,1]),pv_6$fit+pv_6$se.fit*1.96,lwd=1,lty=3)
lines((pd_6[,1]),pv_6$fit-pv_6$se.fit *1.96,lwd=1,lty=3)
text(1,1000,"(a)",cex=1.5)

pd_7 <- data.frame(Week=seq(1,10,length=30),canola=rep("NAM-13",30))
pv_7 <- predict(gm6,newdata=pd_7,type="link",se=TRUE)
plot(data.root$Week,data.root$Total_RL,type="p",xlab="",ylab="",cex.axis=1.5,cex.lab=1.5, lwd=1)
with(data.root[canola=="NAM-13",],points(Week,subset(Total_RL, canola=="NAM-13")))
lines((pd_7[,1]),pv_7$fit,lwd=3,col="blue")
lines((pd_7[,1]),pv_7$fit+pv_7$se.fit*1.96,lwd=1,lty=3)
lines((pd_7[,1]),pv_7$fit-pv_7$se.fit *1.96,lwd=1,lty=3)
text(1,1000,"(b)",cex=1.5)

pd_8 <- data.frame(Week=seq(1,10,length=30),canola=rep("NAM-14",30))
pv_8 <- predict(gm6,newdata=pd_8,type="link",se=TRUE)
plot(data.root$Week,data.root$Total_RL,type="p",xlab="Sampling Week",ylab="Total Root Length (cm)",cex.axis=1.5,cex.lab=1.5, lwd=1)
with(data.root[canola=="NAM-14",],points(Week,subset(Total_RL, canola=="NAM-14")))
lines((pd_8[,1]),pv_8$fit,lwd=3,col="blue")
lines((pd_8[,1]),pv_8$fit+pv_8$se.fit*1.96,lwd=1,lty=3)
lines((pd_8[,1]),pv_8$fit-pv_8$se.fit *1.96,lwd=1,lty=3)
text(1,1000,"(c)",cex=1.5)

pd_9 <- data.frame(Week=seq(1,10,length=30),canola=rep("NAM-17",30))
pv_9 <- predict(gm6,newdata=pd_9,type="link",se=TRUE)
plot(data.root$Week,data.root$Total_RL,type="p",xlab="Sampling Week",ylab="",cex.axis=1.5,cex.lab=1.5, lwd=1)
with(data.root[canola=="NAM-17",],points(Week,subset(Total_RL, canola=="NAM-17")))
lines((pd_9[,1]),pv_9$fit,lwd=3,col="blue")
lines((pd_9[,1]),pv_9$fit+pv_9$se.fit*1.96,lwd=1,lty=3)
lines((pd_9[,1]),pv_9$fit-pv_9$se.fit *1.96,lwd=1,lty=3)
text(1,1000,"(d)",cex=1.5)

par(mfrow = c(2,2),
    oma = c(5,3,2,0),
    mar = c(4,4,2,3))

pd_10 <- data.frame(Week=seq(1,10,length=30),canola=rep("NAM-23",30))
pv_10 <- predict(gm6,newdata=pd_10,type="link",se=TRUE)
plot(data.root$Week,data.root$Total_RL,type="p",xlab="",ylab="Total Root Length (cm)",cex.axis=1.5,cex.lab=1.5, lwd=1)
with(data.root[canola=="NAM-23",],points(Week,subset(Total_RL, canola=="NAM-23")))
lines((pd_10[,1]),pv_10$fit,lwd=3,col="blue")
lines((pd_10[,1]),pv_10$fit+pv_10$se.fit*1.96,lwd=1,lty=3)
lines((pd_10[,1]),pv_10$fit-pv_10$se.fit *1.96,lwd=1,lty=3)
text(1,1000,"(e)",cex=1.5)

pd_11 <- data.frame(Week=seq(1,10,length=30),canola=rep("NAM-30",30))
pv_11 <- predict(gm6,newdata=pd_11,type="link",se=TRUE)
plot(data.root$Week,data.root$Total_RL,type="p",xlab="",ylab="",cex.axis=1.5,cex.lab=1.5, lwd=1)
with(data.root[canola=="NAM-30",],points(Week,subset(Total_RL, canola=="NAM-30")))
lines((pd_11[,1]),pv_11$fit,lwd=3,col="blue")
lines((pd_11[,1]),pv_11$fit+pv_11$se.fit*1.96,lwd=1,lty=3)
lines((pd_11[,1]),pv_11$fit-pv_11$se.fit *1.96,lwd=1,lty=3)
text(1,1000,"(f)",cex=1.5)

pd_12 <- data.frame(Week=seq(1,10,length=30),canola=rep("NAM-32",30))
pv_12 <- predict(gm6,newdata=pd_12,type="link",se=TRUE)
plot(data.root$Week,data.root$Total_RL,type="p",xlab="Sampling Week",ylab="Total Root Length (cm)",cex.axis=1.5,cex.lab=1.5, lwd=1)
with(data.root[canola=="NAM-32",],points(Week,subset(Total_RL, canola=="NAM-32")))
lines((pd_12[,1]),pv_12$fit,lwd=3,col="blue")
lines((pd_12[,1]),pv_12$fit+pv_12$se.fit*1.96,lwd=1,lty=3)
lines((pd_12[,1]),pv_12$fit-pv_12$se.fit *1.96,lwd=1,lty=3)
text(1,1000,"(g)",cex=1.5)

pd_13 <- data.frame(Week=seq(1,10,length=30),canola=rep("NAM-37",30))
pv_13 <- predict(gm6,newdata=pd_13,type="link",se=TRUE)
plot(data.root$Week,data.root$Total_RL,type="p",xlab="Sampling Week",ylab="",cex.axis=1.5,cex.lab=1.5, lwd=1)
with(data.root[canola=="NAM-37",],points(Week,subset(Total_RL, canola=="NAM-37")))
lines((pd_13[,1]),pv_13$fit,lwd=3,col="blue")
lines((pd_13[,1]),pv_13$fit+pv_13$se.fit*1.96,lwd=1,lty=3)
lines((pd_13[,1]),pv_13$fit-pv_13$se.fit *1.96,lwd=1,lty=3)
text(1,1000,"(h)",cex=1.5)

par(mfrow = c(2,2),
    oma = c(5,3,2,0),
    mar = c(4,4,2,3))

pd_14 <- data.frame(Week=seq(1,10,length=30),canola=rep("NAM-43",30))
pv_14 <- predict(gm6,newdata=pd_14,type="link",se=TRUE)
plot(data.root$Week,data.root$Total_RL,type="p",xlab="",ylab="Total Root Length (cm)",cex.axis=1.5,cex.lab=1.5, lwd=1)
with(data.root[canola=="NAM-43",],points(Week,subset(Total_RL, canola=="NAM-43")))
lines((pd_14[,1]),pv_14$fit,lwd=3,col="blue")
lines((pd_14[,1]),pv_14$fit+pv_14$se.fit*1.96,lwd=1,lty=3)
lines((pd_14[,1]),pv_14$fit-pv_14$se.fit *1.96,lwd=1,lty=3)
text(1,1000,"(i)",cex=1.5)


pd_15 <- data.frame(Week=seq(1,10,length=30),canola=rep("NAM-46",30))
pv_15 <- predict(gm6,newdata=pd_15,type="link",se=TRUE)
plot(data.root$Week,data.root$Total_RL,type="p",xlab="",ylab="",cex.axis=1.5,cex.lab=1.5, lwd=1)
with(data.root[canola=="NAM-46",],points(Week,subset(Total_RL, canola=="NAM-46")))
lines((pd_15[,1]),pv_15$fit,lwd=3,col="blue")
lines((pd_15[,1]),pv_15$fit+pv_15$se.fit*1.96,lwd=1,lty=3)
lines((pd_15[,1]),pv_15$fit-pv_15$se.fit *1.96,lwd=1,lty=3)
text(1,1000,"(j)",cex=1.5)


pd_16 <- data.frame(Week=seq(1,10,length=30),canola=rep("NAM-48",30))
pv_16 <- predict(gm6,newdata=pd_16,type="link",se=TRUE)
plot(data.root$Week,data.root$Total_RL,type="p",xlab="Sampling Week",ylab="Total Root Length (cm)",cex.axis=1.5,cex.lab=1.5, lwd=1)
with(data.root[canola=="NAM-48",],points(Week,subset(Total_RL, canola=="NAM-48")))
lines((pd_16[,1]),pv_16$fit,lwd=3,col="blue")
lines((pd_16[,1]),pv_16$fit+pv_16$se.fit*1.96,lwd=1,lty=3)
lines((pd_16[,1]),pv_16$fit-pv_16$se.fit *1.96,lwd=1,lty=3)
text(1,1000,"(k)",cex=1.5)


pd_17 <- data.frame(Week=seq(1,10,length=30),canola=rep("NAM-5",30))
pv_17 <- predict(gm6,newdata=pd_17,type="link",se=TRUE)
plot(data.root$Week,data.root$Total_RL,type="p",xlab="Sampling Week",ylab="",cex.axis=1.5,cex.lab=1.5, lwd=1)
with(data.root[canola=="NAM-5",],points(Week,subset(Total_RL, canola=="NAM-5")))
lines((pd_17[,1]),pv_17$fit,lwd=3,col="blue")
lines((pd_17[,1]),pv_17$fit+pv_17$se.fit*1.96,lwd=1,lty=3)
lines((pd_17[,1]),pv_17$fit-pv_17$se.fit *1.96,lwd=1,lty=3)
text(1,1000,"(l)",cex=1.5)

par(mfrow = c(2,2),
    oma = c(5,3,2,0),
    mar = c(4,4,2,3))

pd_18 <- data.frame(Week=seq(1,10,length=30),canola=rep("NAM-72",30))
pv_18 <- predict(gm6,newdata=pd_18,type="link",se=TRUE)
plot(data.root$Week,data.root$Total_RL,type="p",xlab="",ylab="Total Root Length (cm)",cex.axis=1.5,cex.lab=1.5, lwd=1)
with(data.root[canola=="NAM-72",],points(Week,subset(Total_RL, canola=="NAM-72")))
lines((pd_18[,1]),pv_18$fit,lwd=3,col="blue")
lines((pd_18[,1]),pv_18$fit+pv_18$se.fit*1.96,lwd=1,lty=3)
lines((pd_18[,1]),pv_18$fit-pv_18$se.fit *1.96,lwd=1,lty=3)
text(1,1000,"(m)",cex=1.5)


pd_19 <- data.frame(Week=seq(1,10,length=30),canola=rep("NAM-76",30))
pv_19 <- predict(gm6,newdata=pd_19,type="link",se=TRUE)
plot(data.root$Week,data.root$Total_RL,type="p",xlab="",ylab="",cex.axis=1.5,cex.lab=1.5, lwd=1)
with(data.root[canola=="NAM-76",],points(Week,subset(Total_RL, canola=="NAM-76")))
lines((pd_19[,1]),pv_19$fit,lwd=3,col="blue")
lines((pd_19[,1]),pv_19$fit+pv_19$se.fit*1.96,lwd=1,lty=3)
lines((pd_19[,1]),pv_19$fit-pv_19$se.fit *1.96,lwd=1,lty=3)
text(1,1000,"(n)",cex=1.5)


pd_20 <- data.frame(Week=seq(1,10,length=30),canola=rep("NAM-79",30))
pv_20 <- predict(gm6,newdata=pd_20,type="link",se=TRUE)
plot(data.root$Week,data.root$Total_RL,type="p",xlab="Sampling Week",ylab="Total Root Length (cm)",cex.axis=1.5,cex.lab=1.5, lwd=1)
with(data.root[canola=="NAM-79",],points(Week,subset(Total_RL, canola=="NAM-79")))
lines((pd_20[,1]),pv_20$fit,lwd=3,col="blue")
lines((pd_20[,1]),pv_20$fit+pv_20$se.fit*1.96,lwd=1,lty=3)
lines((pd_20[,1]),pv_20$fit-pv_20$se.fit *1.96,lwd=1,lty=3)
text(1,1000,"(o)",cex=1.5)


pd_21 <- data.frame(Week=seq(1,10,length=30),canola=rep("NAM-94",30))
pv_21 <- predict(gm6,newdata=pd_21,type="link",se=TRUE)
plot(data.root$Week,data.root$Total_RL,type="p",xlab="Sampling Week",ylab="",cex.axis=1.5,cex.lab=1.5, lwd=1)
with(data.root[canola=="NAM-94",],points(Week,subset(Total_RL, canola=="NAM-94")))
lines((pd_21[,1]),pv_21$fit,lwd=3,col="blue")
lines((pd_21[,1]),pv_21$fit+pv_21$se.fit*1.96,lwd=1,lty=3)
lines((pd_21[,1]),pv_21$fit-pv_21$se.fit *1.96,lwd=1,lty=3)
text(1,1000,"(p)",cex=1.5)

#### summary statistics by canola line and growth stage

library(psych)

RL.summary= describeBy(data.root$Total_RL,list(data.root$Canola.Lines,data.root$Growth.Stage),mat=TRUE)
RL.summary= as.data.frame(RL.summary)
write.csv(RL.summary, 'RL.summary.csv')


########

### Aditional Analysis done based on reviewers comments May 2021
## Correlation test between microbial metrics and root trait metrics 
# Taxonomy glomed at Genus level
# Do mantel test using the whole data set, Vegetative stage, flowering stage, maturity stage datasets separately. 


## Load phyloseq object


## Additional analysis for Manuscript II


library(ggplot2)
library(phyloseq)
library(ape)
library(biomformat)
library(microbiome)
library(ggpubr)
library(knitr)
library(dplyr)
library(ampvis2)
library(radiant) ## this is for rownames_to_column function
library(tibble)## for as.table function
library(microbiomeSeq)
library(metagMisc)

# Set Working directory
setwd("D:/out_puts/phyloseq/with_green_genes_taxa/Additional_analysis_manuscript_II_2021")

# load the phyloseq object containing the  2016 rhizosphere bacterial dataset

load("D:/out_puts/phyloseq/with_green_genes_taxa/2016.physeq.updated.line.names.new.variable.RData")

# Look summary of the phyloseq object

physeq.updated.line.names.new.variable

#### subset dominat genera based on both frequence and mean relative abundance plus previously identified PGPB


physeq.top.10.genera.root.trait =  subset_taxa(physeq.updated.line.names.new.variable, Genus=="Pseudomonas"| Genus=="Flavobacterium"| Genus=="Lysobacter"| Genus=="Rhodoplanes"| Genus=="Pedobacter"| Genus=="Vibrio"| Genus=="Stenotrophomonas"| Genus=="Arthrobacter"| Genus=="Gluconacetobacter"| Genus=="Skermanella"| Genus=="Bradyrhizobium"| Genus=="Paenibacillus"| Genus=="Pantoea"| Genus=="Rhizobium"| Genus=="Bacillus"| Genus=="Leifsonia"| Genus=="Microbacterium"| Genus=="Rhodococcus"| Genus=="Xanthomonas"| Genus=="Variovorax"| Genus=="Agrobacterium"| Genus=="Phyllobacterium")

## remove zero sum taxa
physeq.top.10.genera.root.trait= prune_taxa(taxa_sums(physeq.top.10.genera.root.trait) > 0, physeq.top.10.genera.root.trait)


## remove zero sum samples
#physeq.top.10.genera.root.trait = prune_samples(sample_sums(physeq.top.10.genera.root.trait)>=1, physeq.top.10.genera.root.trait)

#physeq <- taxa_level(physeq.top.10.genera.root.trait, "Genus") ### tax_glom at geneus level with genus name indicated

physeq_genus<- tax_glom(physeq.top.10.genera.root.trait, taxrank="Genus")

## remove zero sum samples
physeq_genus = prune_samples(sample_sums(physeq_genus)>=1, physeq_genus)

# Change the phyloseq to CSV

physeq_melt <- psmelt(physeq_genus)
write.csv(physeq_melt, 'physeq_melt.csv')

## open the physeq_melt.csv file and prepare it for mantel test

# load the csv file with subset genera and metadata
Most_abundant_pgpb_subset_mantel_correlation_analysis

## run mantel test on bacterial matrix vs different root traits and different growth stages

# sunset the data for preparing the distace matrix

library(vegan)
# substet variable of interest: using all data set
abund = Most_abundant_pgpb_subset_mantel_correlation_analysis[,26:ncol(Most_abundant_pgpb_subset_mantel_correlation_analysis)] # bacterial genera begins starting column 26

root_length_total = Most_abundant_pgpb_subset_mantel_correlation_analysis$Total_RL

fine_root_length_total = Most_abundant_pgpb_subset_mantel_correlation_analysis$TFRL2mm
coarse_root_length_total = Most_abundant_pgpb_subset_mantel_correlation_analysis$CoarseRL
soil_moisture = Most_abundant_pgpb_subset_mantel_correlation_analysis$Soil.moisture

FR = Most_abundant_pgpb_subset_mantel_correlation_analysis$FRL

enviro_variables = Most_abundant_pgpb_subset_mantel_correlation_analysis[,20:25]
extra_fine_root_length = Most_abundant_pgpb_subset_mantel_correlation_analysis$ExtraFRL
rbiomass = Most_abundant_pgpb_subset_mantel_correlation_analysis$Root.Biomass

##
#abundance data frame - bray curtis dissimilarity
dist.abund = vegdist(abund, method = "bray")

#root length total - euclidean distance
dist.root_length_total = dist(root_length_total, method = "euclidean")

#fine root length total - euclidean distance

dist.fine_root_length_total = dist(fine_root_length_total, method = "euclidean")

#fine root length  - euclidean distance

dist.FRL = dist(FR, method = "euclidean")

# coarse root length- euclidean
dist.coarse_root_length_total = dist(coarse_root_length_total, method = "euclidean")

# root biomass
dist.rbiomass= dist(rbiomass, method = "euclidean")


#soil moisture data frame - euclidean 
dist.soil_moisture = dist(soil_moisture, method = "euclidean")


# all environmental variable data frame - euclidean 
dist.enviro_variables = dist(scale(enviro_variables), method = "euclidean")

# Extra fine root length  - euclidean distance

dist.extra_fine_root_length = dist(extra_fine_root_length, method = "euclidean")

#abundance vs biomass
abund_rbiomass = mantel(dist.abund, dist.rbiomass, method = "spearman", permutations = 1000, na.rm = TRUE)
abund_rbiomass

#abundance vs total root length 
abund_TRL = mantel(dist.abund, dist.root_length_total, method = "spearman", permutations = 1000, na.rm = TRUE)
abund_TRL

# abundance vs total fine root length
abund_TFRL = mantel(dist.abund, dist.fine_root_length_total, method = "spearman", permutations = 1000, na.rm = TRUE)
abund_TFRL

# abundance vs coarse  root length
abund_TCRL = mantel(dist.abund, dist.coarse_root_length_total, method = "spearman", permutations = 1000, na.rm = TRUE)
abund_TCRL

#abundane vs soil moisuture

abund_dist.soil_moisture = mantel(dist.abund, dist.soil_moisture, method = "spearman", permutations = 1000, na.rm = TRUE)
abund_dist.soil_moisture

# abundance vs all environmental variables

abund_dist.enviro_variables = mantel(dist.abund, dist.enviro_variables, method = "spearman", permutations = 1000, na.rm = TRUE)
abund_dist.enviro_variables

# abundance vs extra fine root length
abund_EFRL = mantel(dist.abund, dist.extra_fine_root_length, method = "spearman", permutations = 1000, na.rm = TRUE)
abund_EFRL

# abundance vs  fine root length
abund_FRL = mantel(dist.abund, dist.FRL, method = "spearman", permutations = 1000, na.rm = TRUE)
abund_FRL

## Subset vegetative growth stage data

Vegetative_data <- subset(Most_abundant_pgpb_subset_mantel_correlation_analysis,Growth.Stage =='Vegetative')

vabund = Vegetative_data[,26:ncol(Vegetative_data)] # bacterial genera begins starting column 26

vroot_length_total = Vegetative_data$Total_RL

vfine_root_length_total = Vegetative_data$TFRL2mm
vcoarse_root_length_total = Vegetative_data$CoarseRL
vFR = Vegetative_data$FRL
vEFR= Vegetative_data$ExtraFRL
vrbiomass = Vegetative_data$Root.Biomass
vsoil_moisture = Vegetative_data$Soil.moisture

venviro_variables = Vegetative_data[,20:25]

#abundance data frame - bray curtis dissimilarity
vdist.abund = vegdist(vabund, method = "bray")

#root biomass - euclidean distance
vdist.vrbiomass = dist(vrbiomass, method = "euclidean")

#root length total - euclidean distance
vdist.root_length_total = dist(vroot_length_total, method = "euclidean")

#fine root length total - euclidean distance

vdist.fine_root_length_total = dist(vfine_root_length_total, method = "euclidean")

#fine root length  - euclidean distance

vdist.FRL = dist(vFR, method = "euclidean")

# Extra fine root length  - euclidean distance

vdist.EFRL = dist(vEFR, method = "euclidean")

# coarse root length- euclidean
vdist.coarse_root_length_total = dist(vcoarse_root_length_total, method = "euclidean")



#soil moisture data frame - euclidean 
vdist.soil_moisture = dist(vsoil_moisture, method = "euclidean")


# all environmental variable data frame - euclidean 
vdist.enviro_variables = dist(venviro_variables, method = "euclidean")


#abundance vs root biomass
vabund_rbiomass = mantel(vdist.abund, vdist.vrbiomass, method = "spearman", permutations = 1000, na.rm = TRUE)
vabund_rbiomass

#abundance vs total root length 
vabund_TRL = mantel(vdist.abund, vdist.root_length_total, method = "spearman", permutations = 1000, na.rm = TRUE)
vabund_TRL

# abundance vs total fine root length
vabund_TFRL = mantel(vdist.abund, vdist.fine_root_length_total, method = "spearman", permutations = 1000, na.rm = TRUE)
vabund_TFRL


# abundance vs  fine root length
vabund_FRL = mantel(vdist.abund, vdist.FRL, method = "spearman", permutations = 1000, na.rm = TRUE)
vabund_FRL

# abundance vs  Extra fine root length
vabund_EFRL = mantel(vdist.abund, vdist.EFRL, method = "spearman", permutations = 1000, na.rm = TRUE)
vabund_EFRL

# abundance vs coarse  root length
vabund_TCRL = mantel(vdist.abund, vdist.coarse_root_length_total, method = "spearman", permutations = 1000, na.rm = TRUE)
vabund_TCRL

#abundane vs soil moisuture

vabund_dist.soil_moisture = mantel(vdist.abund, vdist.soil_moisture, method = "spearman", permutations = 1000, na.rm = TRUE)
vabund_dist.soil_moisture

# abundance vs all environmental variables

vabund_dist.enviro_variables = mantel(vdist.abund, vdist.enviro_variables, method = "spearman", permutations = 1000, na.rm = TRUE)
vabund_dist.enviro_variables

########

## Subset Flowering growth stage data

Flowering_data <- subset(Most_abundant_pgpb_subset_mantel_correlation_analysis,Growth.Stage =='Flowering')

fabund = Flowering_data[,26:ncol(Flowering_data)] # bacterial genera begins starting column 26

froot_length_total = Flowering_data$Total_RL

ffine_root_length_total = Flowering_data$TFRL2mm
fcoarse_root_length_total = Flowering_data$CoarseRL
fFRL = Flowering_data$FRL
fEFRL = Flowering_data$ExtraFRL
frbiomass = Flowering_data$Root.Biomass

fsoil_moisture = Flowering_data$Soil.moisture

fenviro_variables = Flowering_data[,20:25]

#abundance data frame - bray curtis dissimilarity
fdist.abund = vegdist(fabund, method = "bray")

#root biomass  - euclidean distance
fdist.rbiomass = dist(frbiomass, method = "euclidean")

#root length total - euclidean distance
fdist.root_length_total = dist(froot_length_total, method = "euclidean")

#fine root length total - euclidean distance

fdist.fine_root_length_total = dist(ffine_root_length_total, method = "euclidean")

#fine root length  - euclidean distance

fdist.FRL = dist(fFRL, method = "euclidean")

# Extra fine root length  - euclidean distance

fdist.EFRL = dist(fEFRL, method = "euclidean")


# coarse root length- euclidean
fdist.coarse_root_length_total = dist(fcoarse_root_length_total, method = "euclidean")



#soil moisture data frame - euclidean 
fdist.soil_moisture = dist(fsoil_moisture, method = "euclidean")


# all environmental variable data frame - euclidean 
fdist.enviro_variables = dist(fenviro_variables, method = "euclidean")

#abundance vs  root biomass 
fabund_rbiomass = mantel(fdist.abund, fdist.rbiomass, method = "spearman", permutations = 1000, na.rm = TRUE)
fabund_rbiomass


#abundance vs total root length 
fabund_TRL = mantel(fdist.abund, fdist.root_length_total, method = "spearman", permutations = 1000, na.rm = TRUE)
fabund_TRL

# abundance vs total fine root length
fabund_TFRL = mantel(fdist.abund, fdist.fine_root_length_total, method = "spearman", permutations = 1000, na.rm = TRUE)
fabund_TFRL


# abundance vs  fine root length
fabund_FRL = mantel(fdist.abund, fdist.FRL, method = "spearman", permutations = 1000, na.rm = TRUE)
fabund_FRL

# abundance vs  Extra fine root length
fabund_EFRL = mantel(fdist.abund, fdist.EFRL, method = "spearman", permutations = 1000, na.rm = TRUE)
fabund_EFRL

# abundance vs coarse  root length
fabund_TCRL = mantel(fdist.abund, fdist.coarse_root_length_total, method = "spearman", permutations = 1000, na.rm = TRUE)
fabund_TCRL

#abundane vs soil moisuture

fabund_dist.soil_moisture = mantel(fdist.abund, fdist.soil_moisture, method = "spearman", permutations = 1000, na.rm = TRUE)
fabund_dist.soil_moisture

# abundance vs all environmental variables

fabund_dist.enviro_variables = mantel(fdist.abund, fdist.enviro_variables, method = "spearman", permutations = 1000, na.rm = TRUE)
fabund_dist.enviro_variables

#### Partial mantel test controlling for sampling week and genotypes
library(plyr)
# revalue the canola line character values to numerical to calculate distance

Most_abundant_pgpb_subset_mantel_correlation_analysis$genotype <- revalue(Most_abundant_pgpb_subset_mantel_correlation_analysis$Canola.Lines, c("NAM-0"="1", "NAM-13"="2", "NAM-14"="3", "NAM-17"="4", "NAM-5"="5", "NAM-23"="6", "NAM-30"="7", "NAM-32"="8", "NAM-37"="9", "NAM-43"="10", "NAM-46"="11", "NAM-48"="12", "NAM-72"="13", "NAM-76"="14", "NAM-79"="15","NAM-94"="16"))

# select growth stage and sampling week
library(tidyverse)

# all data 
weeek_canola = Most_abundant_pgpb_subset_mantel_correlation_analysis %>% select(genotype, Week)

dis_week_canola = dist(weeek_canola, method = "euclidean")


#abundance vs biomass controlling for sampling week and canola genotype

abund_rbiomass_wc = mantel.partial(dist.abund, dist.rbiomass,dis_week_canola, method = "spearman", permutations = 1000, na.rm = TRUE)
abund_rbiomass_wc

# abundance vs coarse  root length controlling for sampling week and canola genotype


abund_TCRL_wc = mantel.partial(dist.abund, dist.coarse_root_length_total,dis_week_canola, method = "pearson", permutations = 1000, na.rm = TRUE)
abund_TCRL_wc


#abundance vs total root length 
abund_TRL_wc = mantel.partial(dist.abund, dist.root_length_total,dis_week_canola, method = "spearman", permutations = 1000, na.rm = TRUE)
abund_TRL_wc

# abundance vs total fine root length
abund_TFRL_wc = mantel.partial(dist.abund, dist.fine_root_length_total,dis_week_canola, method = "spearman", permutations = 1000, na.rm = TRUE)
abund_TFRL_wc


# abundance vs  fine root length
abund_FRL_wc = mantel.partial(dist.abund, dist.FRL,dis_week_canola, method = "spearman", permutations = 1000, na.rm = TRUE)
abund_FRL_wc

# abundance vs  Extra fine root length
abund_EFRL_wc = mantel.partial(dist.abund, dist.extra_fine_root_length,dis_week_canola, method = "spearman", permutations = 1000, na.rm = TRUE)
abund_EFRL_wc

#abundane vs soil moisuture

abund_dist.soil_moisture_wc = mantel.partial(dist.abund, dist.soil_moisture, dis_week_canola, method = "spearman", permutations = 1000, na.rm = TRUE)
abund_dist.soil_moisture_wc

# abundance vs all environmental variables

abund_dist.enviro_variables_wc = mantel.partial(dist.abund, dist.enviro_variables,dis_week_canola, method = "spearman", permutations = 1000, na.rm = TRUE)
abund_dist.enviro_variables_wc

## Vegetative stage controlling for sampling week and genotype

Vegetative_data <- subset(Most_abundant_pgpb_subset_mantel_correlation_analysis,Growth.Stage =='Vegetative')

vweeek_canola = Vegetative_data %>% select(genotype, Week)
vdis_week_canola = dist(vweeek_canola, method = "euclidean")

#

#abundance vs root biomass
vabund_rbiomass_wc = mantel.partial(vdist.abund, vdist.vrbiomass,vdis_week_canola, method = "spearman", permutations = 1000, na.rm = TRUE)
vabund_rbiomass_wc
#abundance vs total root length 
vabund_TRL_wc = mantel.partial(vdist.abund, vdist.root_length_total,vdis_week_canola, method = "spearman", permutations = 1000, na.rm = TRUE)
vabund_TRL_wc

# abundance vs total fine root length
vabund_TFRL_wc = mantel.partial(vdist.abund, vdist.fine_root_length_total,vdis_week_canola, method = "spearman", permutations = 1000, na.rm = TRUE)
vabund_TFRL_wc


# abundance vs  fine root length
vabund_FRL_wc = mantel.partial(vdist.abund, vdist.FRL,vdis_week_canola, method = "spearman", permutations = 1000, na.rm = TRUE)
vabund_FRL_wc

# abundance vs  Extra fine root length
vabund_EFRL_wc = mantel.partial(vdist.abund, vdist.EFRL,vdis_week_canola, method = "spearman", permutations = 1000, na.rm = TRUE)
vabund_EFRL_wc



# abundance vs coarse  root length
vabund_TCRL_wc = mantel.partial(vdist.abund, vdist.coarse_root_length_total,vdis_week_canola, method = "spearman", permutations = 1000, na.rm = TRUE)
vabund_TCRL_wc

#abundane vs soil moisuture

vabund_dist.soil_moisture_wc = mantel.partial(vdist.abund, vdist.soil_moisture, vdis_week_canola, method = "spearman", permutations = 1000, na.rm = TRUE)
vabund_dist.soil_moisture_wc

# abundance vs all environmental variables

vabund_dist.enviro_variables_wc = mantel.partial(vdist.abund, vdist.enviro_variables,vdis_week_canola, method = "spearman", permutations = 1000, na.rm = TRUE)
vabund_dist.enviro_variables_wc

#######F

## Flowering stage controlling for sampling week and genotype

Flowering_data <- subset(Most_abundant_pgpb_subset_mantel_correlation_analysis,Growth.Stage =='Flowering')

fweeek_canola = Flowering_data %>% select(genotype, Week)
fdis_week_canola = dist(fweeek_canola, method = "euclidean")

#
#abundance vs  root biomass 
fabund_rbiomass_wc = mantel.partial(fdist.abund, fdist.rbiomass,fdis_week_canola, method = "spearman", permutations = 1000, na.rm = TRUE)
fabund_rbiomass_wc

#abundance vs total root length 
fabund_TRL_wc = mantel.partial(fdist.abund, fdist.root_length_total,fdis_week_canola, method = "spearman", permutations = 1000, na.rm = TRUE)
fabund_TRL_wc

# abundance vs total fine root length
fabund_TFRL_wc = mantel.partial(fdist.abund, fdist.fine_root_length_total,fdis_week_canola, method = "spearman", permutations = 1000, na.rm = TRUE)
fabund_TFRL_wc


# abundance vs  fine root length
fabund_FRL_wc = mantel.partial(fdist.abund, fdist.FRL,fdis_week_canola, method = "spearman", permutations = 1000, na.rm = TRUE)

fabund_FRL_wc

# abundance vs  Extra fine root length
fabund_EFRL_wc = mantel.partial(fdist.abund, fdist.EFRL,fdis_week_canola, method = "spearman", permutations = 1000, na.rm = TRUE)

fabund_EFRL_wc


# abundance vs coarse  root length
fabund_TCRL_wc = mantel.partial(fdist.abund, fdist.coarse_root_length_total,fdis_week_canola, method = "spearman", permutations = 1000, na.rm = TRUE)
fabund_TCRL_wc

#abundane vs soil moisuture

fabund_dist.soil_moisture_wc = mantel.partial(fdist.abund, fdist.soil_moisture, fdis_week_canola, method = "spearman", permutations = 1000, na.rm = TRUE)
fabund_dist.soil_moisture_wc

# abundance vs all environmental variables

fabund_dist.enviro_variables_wc = mantel.partial(fdist.abund, fdist.enviro_variables,fdis_week_canola, method = "spearman", permutations = 1000, na.rm = TRUE)
fabund_dist.enviro_variables_wc

#### correlation between microbial composition and canola yield per plot 

# using all dataset
canola_yield_plot = Most_abundant_pgpb_subset_mantel_correlation_analysis$Yield.per.plot

vcanola_yield_plot = Vegetative_data$Yield.per.plot
fcanola_yield_plot = Flowering_data$Yield.per.plot


dis_canola_yield_plot= dist(canola_yield_plot, method = "euclidean")
vdis_vcanola_yield_plot= dist(vcanola_yield_plot, method = "euclidean")
fdis_fcanola_yield_plot= dist(fcanola_yield_plot, method = "euclidean")


#abundance vs yield

abund_yield = mantel(dist.abund, dis_canola_yield_plot, method = "spearman", permutations = 1000, na.rm = TRUE)
abund_yield

vabund_yield = mantel(vdist.abund, vdis_vcanola_yield_plot, method = "spearman", permutations = 1000, na.rm = TRUE)
vabund_yield

fabund_yield = mantel(fdist.abund, fdis_fcanola_yield_plot, method = "spearman", permutations = 1000, na.rm = TRUE)
fabund_yield

#root legth vs yield

RL_yield = mantel(dist.root_length_total, dis_canola_yield_plot, method = "spearman", permutations = 1000, na.rm = TRUE)
RL_yield

vRL_yield = mantel(vdist.root_length_total, vdis_vcanola_yield_plot, method = "spearman", permutations = 1000, na.rm = TRUE)
vRL_yield

fRL_yield = mantel(fdist.root_length_total, fdis_fcanola_yield_plot, method = "spearman", permutations = 1000, na.rm = TRUE)
fRL_yield


# Fine root legth vs yield

RL_yield = mantel(dist.fine_root_length_total, dis_canola_yield_plot, method = "spearman", permutations = 1000, na.rm = TRUE)
RL_yield

vRL_yield = mantel(vdist.fine_root_length_total, vdis_vcanola_yield_plot, method = "spearman", permutations = 1000, na.rm = TRUE)
vRL_yield

fRL_yield = mantel(fdist.fine_root_length_total, fdis_fcanola_yield_plot, method = "spearman", permutations = 1000, na.rm = TRUE)
fRL_yield




####

#### Correlation between Root length at week 1 and 2 with yield

# subset week 1, 2 and 1 and 2 data

Week_1_2_data <- subset(Most_abundant_pgpb_subset_mantel_correlation_analysis,Week =='1'|Week=='2')
Week_1_data <- subset(Most_abundant_pgpb_subset_mantel_correlation_analysis,Week =='1')
Week_2_data <- subset(Most_abundant_pgpb_subset_mantel_correlation_analysis,Week =='2')

#
# substet taxa variable of interest: w1
abund_w1 = Week_1_data[,26:46] # bacterial genera begins starting column 26
abund_w2 = Week_2_data[,26:46] # bacterial genera begins starting column 26
abund_w1_2 = Week_1_2_data[,26:46] # bacterial genera begins starting column 26

#abundance data frame - bray curtis dissimilarity
w1_dist.abund = vegdist(abund_w1, method = "bray")
w2_dist.abund = vegdist(abund_w2, method = "bray")
w1_2_dist.abund = vegdist(abund_w1_2, method = "bray")


w1_canola_yield_plot = Week_1_data$Yield.per.plot
w2_canola_yield_plot = Week_2_data$Yield.per.plot
w1_2_canola_yield_plot = Week_1_2_data$Yield.per.plot


w1_dis_canola_yield_plot= dist(w1_canola_yield_plot, method = "euclidean")
w2_dis_canola_yield_plot= dist(w2_canola_yield_plot, method = "euclidean")
w1_2_dis_canola_yield_plot= dist(w1_2_canola_yield_plot, method = "euclidean")

##
w1_root_length_total = Week_1_data$Total_RL
w2_root_length_total = Week_2_data$Total_RL
w1_2_root_length_total = Week_1_2_data$Total_RL

w1_TFRL = Week_1_data$TFRL2mm
w2_TFRL = Week_2_data$TFRL2mm
w1_2_TFRL = Week_1_2_data$TFRL2mm

w1_FRL = Week_1_data$FRL
w2_FRL = Week_2_data$FRL
w1_2_FRL = Week_1_2_data$FRL

w1_EFRL = Week_1_data$ExtraFRL
w2_EFRL = Week_2_data$ExtraFRL
w1_2_EFRL = Week_1_2_data$ExtraFRL

w1_CRL = Week_1_data$CoarseRL
w2_CRL = Week_2_data$CoarseRL
w1_2_CRL = Week_1_2_data$CoarseRL

#
w1_dis_root_length_total= dist(w1_root_length_total, method = "euclidean")
w2_dis_root_length_total= dist(w2_root_length_total, method = "euclidean")
w1_2_dis_root_length_total= dist(w1_2_root_length_total, method = "euclidean")

w1_dis_TFR= dist(w1_TFRL, method = "euclidean")
w2_dis_TFR= dist(w2_TFRL, method = "euclidean")
w1_2_TFR= dist(w1_2_TFRL, method = "euclidean")

w1_dis_FR= dist(w1_FRL, method = "euclidean")
w2_dis_FR= dist(w2_FRL, method = "euclidean")
w1_2_FR= dist(w1_2_FRL, method = "euclidean")

w1_dis_EFR= dist(w1_EFRL, method = "euclidean")
w2_dis_EFR= dist(w2_EFRL, method = "euclidean")
w1_2_EFR= dist(w1_2_EFRL, method = "euclidean")


w1_dis_CR= dist(w1_CRL, method = "euclidean")
w2_dis_CR= dist(w2_CRL, method = "euclidean")
w1_2_CR= dist(w1_2_CRL, method = "euclidean")

# Total root legth vs yield: not significant

w1_RL_yield = mantel(w1_dis_root_length_total, w1_dis_canola_yield_plot, method = "spearman", permutations = 1000, na.rm = TRUE)
w1_RL_yield

w2_RL_yield = mantel(w2_dis_root_length_total, w2_dis_canola_yield_plot, method = "spearman", permutations = 1000, na.rm = TRUE)
w2_RL_yield

w1_2_RL_yield = mantel(w1_2_dis_root_length_total, w1_2_dis_canola_yield_plot, method = "spearman", permutations = 1000, na.rm = TRUE)
w1_2_RL_yield

# Total fine root legth vs yield : NOt significant

w1_FRL_yield = mantel(w1_dis_TFR, w1_dis_canola_yield_plot, method = "spearman", permutations = 1000, na.rm = TRUE)
w1_FRL_yield

w2_FRL_yield = mantel(w2_dis_TFR, w2_dis_canola_yield_plot, method = "spearman", permutations = 1000, na.rm = TRUE)
w2_FRL_yield

w1_2_FRL_yield = mantel(w1_2_TFR, w1_2_dis_canola_yield_plot, method = "spearman", permutations = 1000, na.rm = TRUE)
w1_2_FRL_yield


# fine root legth vs yield :  significant

w1_FR_yield = mantel(w1_dis_FR, w1_dis_canola_yield_plot, method = "spearman", permutations = 1000, na.rm = TRUE)
w1_FR_yield  ## significant correlation at week1

# here lets control for genotype difference to see how much fine root length at w1 is correlated with canola yield

w1_canola = Week_1_data %>% select(genotype) # subset genotype
w1_dis_canola = dist(w1_canola, method = "euclidean")

# No run mantel test controlling for genotype difference

# fine root legth vs yield :  significant

w1_FR_yield_genotype = mantel.partial(w1_dis_FR, w1_dis_canola_yield_plot,w1_dis_canola, method = "spearman", permutations = 1000, na.rm = TRUE)
w1_FR_yield_genotype 

#

w2_FR_yield = mantel(w2_dis_FR, w2_dis_canola_yield_plot, method = "spearman", permutations = 1000, na.rm = TRUE)
w2_FR_yield

w1_2_FR_yield = mantel(w1_2_FR, w1_2_dis_canola_yield_plot, method = "spearman", permutations = 1000, na.rm = TRUE)
w1_2_FR_yield


# Extra fine root legth vs yield : NOt significant

w1_EFR_yield = mantel(w1_dis_EFR, w1_dis_canola_yield_plot, method = "spearman", permutations = 1000, na.rm = TRUE)
w1_EFR_yield  

w2_EFR_yield = mantel(w2_dis_EFR, w2_dis_canola_yield_plot, method = "spearman", permutations = 1000, na.rm = TRUE)
w2_EFR_yield

w1_2_EFR_yield = mantel(w1_2_EFR, w1_2_dis_canola_yield_plot, method = "spearman", permutations = 1000, na.rm = TRUE)
w1_2_EFR_yield

# Coarse root legth vs yield : NOt significant

w1_CR_yield = mantel(w1_dis_CR, w1_dis_canola_yield_plot, method = "spearman", permutations = 1000, na.rm = TRUE)
w1_CR_yield  

w2_CR_yield = mantel(w2_dis_CR, w2_dis_canola_yield_plot, method = "spearman", permutations = 1000, na.rm = TRUE)
w2_CR_yield

w1_2_CR_yield = mantel(w1_2_CR, w1_2_dis_canola_yield_plot, method = "spearman", permutations = 1000, na.rm = TRUE)
w1_2_CR_yield


##### mantel test on week 1 data where FRL showed significant correlation with yield 

w1_dist.abund
w1_dis_FR
w1_dis_canola

w1_abund_FR_genotype = mantel.partial(w1_dist.abund, w1_dis_FR,w1_dis_canola, method = "spearman", permutations = 1000, na.rm = TRUE)
w1_abund_FR_genotype



## All data Controlling for environmental variables 


abund_rbiomass_ev = mantel.partial(dist.abund, dist.rbiomass,dist.enviro_variables, method = "spearman", permutations = 1000, na.rm = TRUE)
abund_rbiomass_ev

# abundance vs coarse  root length controlling for environmental variables

abund_TCRL_ev = mantel.partial(dist.abund, dist.coarse_root_length_total,dist.enviro_variables, method = "pearson", permutations = 1000, na.rm = TRUE)
abund_TCRL_ev


#abundance vs total root length 
abund_TRL_ev = mantel.partial(dist.abund, dist.root_length_total,dist.enviro_variables, method = "spearman", permutations = 1000, na.rm = TRUE)
abund_TRL_ev

# abundance vs total fine root length
abund_TFRL_ev = mantel.partial(dist.abund, dist.fine_root_length_total,dist.enviro_variables, method = "spearman", permutations = 1000, na.rm = TRUE)
abund_TFRL_ev


# abundance vs  fine root length
abund_FRL_ev = mantel.partial(dist.abund, dist.FRL,dist.enviro_variables, method = "spearman", permutations = 1000, na.rm = TRUE)
abund_FRL_ev

# abundance vs  Extra fine root length
abund_EFRL_ev = mantel.partial(dist.abund, dist.extra_fine_root_length,dist.enviro_variables, method = "spearman", permutations = 1000, na.rm = TRUE)
abund_EFRL_ev

# envi vs coarse  root length
fTCRL_ev = mantel(dist.root_length_total,dist.enviro_variables, method = "spearman", permutations = 1000, na.rm = TRUE)

fTCRL_ev

#######
## Vegetative stage data: controlling for environmental variables


#abundance vs root biomass
vabund_rbiomass_ev = mantel.partial(vdist.abund, vdist.vrbiomass,vdist.enviro_variables, method = "spearman", permutations = 1000, na.rm = TRUE)
vabund_rbiomass_ev

#abundance vs total root length 
vabund_TRL_ev = mantel.partial(vdist.abund, vdist.root_length_total,vdist.enviro_variables, method = "spearman", permutations = 1000, na.rm = TRUE)
vabund_TRL_ev

# abundance vs total fine root length
vabund_TFRL_ev = mantel.partial(vdist.abund, vdist.fine_root_length_total,vdist.enviro_variables, method = "spearman", permutations = 1000, na.rm = TRUE)
vabund_TFRL_ev


# abundance vs  fine root length
vabund_FRL_ev = mantel.partial(vdist.abund, vdist.FRL,vdist.enviro_variables, method = "spearman", permutations = 1000, na.rm = TRUE)
vabund_FRL_ev

# abundance vs  Extra fine root length
vabund_EFRL_ev = mantel.partial(vdist.abund, vdist.EFRL,vdist.enviro_variables, method = "spearman", permutations = 1000, na.rm = TRUE)
vabund_EFRL_ev


# abundance vs coarse  root length
vabund_TCRL_ev = mantel.partial(vdist.abund, vdist.coarse_root_length_total,vdist.enviro_variables, method = "spearman", permutations = 1000, na.rm = TRUE)
vabund_TCRL_ev

# envi vs coarse  root length
fTCRL_ev = mantel(vdist.root_length_total,vdist.enviro_variables, method = "spearman", permutations = 1000, na.rm = TRUE)

fTCRL_ev
##################

# Flowering stage: controlling for environmental variables
#abundance vs  root biomass 

fabund_rbiomass_ev = mantel.partial(fdist.abund, fdist.rbiomass,fdist.enviro_variables, method = "spearman", permutations = 1000, na.rm = TRUE)
fabund_rbiomass_ev

#abundance vs total root length 
fabund_TRL_ev = mantel.partial(fdist.abund, fdist.root_length_total,fdist.enviro_variables, method = "spearman", permutations = 1000, na.rm = TRUE)

fabund_TRL_ev

# abundance vs total fine root length
fabund_TFRL_ev = mantel.partial(fdist.abund, fdist.fine_root_length_total,fdist.enviro_variables, method = "spearman", permutations = 1000, na.rm = TRUE)

fabund_TFRL_ev


# abundance vs  fine root length
fabund_FRL_ev = mantel.partial(fdist.abund, fdist.FRL,fdist.enviro_variables, method = "spearman", permutations = 1000, na.rm = TRUE)

fabund_FRL_ev

# abundance vs  Extra fine root length
fabund_EFRL_ev = mantel.partial(fdist.abund, fdist.EFRL,fdist.enviro_variables, method = "spearman", permutations = 1000, na.rm = TRUE)

fabund_EFRL_ev


# abundance vs coarse  root length
fabund_TCRL_ev = mantel.partial(fdist.abund, fdist.coarse_root_length_total,fdist.enviro_variables, method = "spearman", permutations = 1000, na.rm = TRUE)

fabund_TCRL_ev



# envi vs coarse  root length
fTCRL_ev = mantel(fdist.root_length_total,fdist.enviro_variables, method = "spearman", permutations = 1000, na.rm = TRUE)

fTCRL_ev

## Additional Analysis
# Linear regression: Bacteria by root traits

library(vegan)

mydata=Root_bacteria

#distribution

hist(log(mydata$Agrobacterium)) # normal dis
hist(log(mydata$Arthrobacter))  # normal
hist(log(mydata$Bacillus))# some zero

hist(log(mydata$Bradyrhizobium))#normal
hist(log(mydata$Flavobacterium))#normal
hist(log(mydata$Gluconacetobacter))#normal
hist(log(mydata$Lysobacter))#normal
hist(log(mydata$Microbacterium))#normal
hist(log(mydata$Paenibacillus))#normal
hist(log(mydata$Pedobacter))#normal
hist(log(mydata$Phyllobacterium))#normal
hist(log(mydata$Pseudomonas))#normal
hist(log(mydata$Rhizobium))#normal
hist(log(mydata$Rhodococcus))#normal

hist(log(mydata$Rhodoplanes))#normal
hist(log(mydata$Skermanella))#normal
hist(log(mydata$Stenotrophomonas))#normal
hist(log(mydata$Variovorax))#normal
hist(log(mydata$Xanthomonas))#normal

#Fit linear models on log transformed bacterial genera abudnance and root traits
  
## Agrobacterium

m1 = lm((log(mydata$Agrobacterium + 1))~log(mydata$Total_RL), data=mydata)
m1_1 = lm((log(mydata$Agrobacterium + 1))~log(mydata$TFRL2mm), data=mydata)
m1_2 = lm((log(mydata$Agrobacterium + 1))~log(mydata$FRL + 1), data=mydata)
m1_3 = lm((log(mydata$Agrobacterium + 1))~log(mydata$ExtraFRL), data=mydata)
m1_4 = lm((log(mydata$Agrobacterium + 1))~log(mydata$CoarseRL + 1), data=mydata)

#Summary
summary(m1)
summary(m1_1)
summary(m1_2)
summary(m1_3)
summary(m1_4)


## Arthrobacter

m2 = lm((log(mydata$Arthrobacter + 1))~log(mydata$Total_RL), data=mydata)
m2_1 = lm((log(mydata$Arthrobacter + 1))~log(mydata$TFRL2mm), data=mydata)
m2_2 = lm((log(mydata$Arthrobacter + 1))~log(mydata$FRL +1), data=mydata)
m2_3 = lm((log(mydata$Arthrobacter + 1))~log(mydata$ExtraFRL), data=mydata)
m2_4 = lm((log(mydata$Arthrobacter + 1))~log(mydata$CoarseRL + 1), data=mydata)

#Summary
summary(m2)
summary(m2_1)
summary(m2_2)
summary(m2_3)
summary(m2_4)

## Bacillus
m3 = lm((log(mydata$Bacillus + 1))~log(mydata$Total_RL), data=mydata)
m3_1 = lm((log(mydata$Bacillus + 1))~log(mydata$TFRL2mm), data=mydata)
m3_2 = lm((log(mydata$Bacillus + 1))~log(mydata$FRL + 1), data=mydata)
m3_3 = lm((log(mydata$Bacillus + 1))~log(mydata$ExtraFRL), data=mydata)
m3_4 = lm((log(mydata$Bacillus + 1))~log(mydata$CoarseRL + 1), data=mydata)

#Summary
summary(m3)
summary(m3_1)
summary(m3_2)
summary(m3_3)
summary(m3_4)

##Bradyrhizobium
m4 = lm((log(mydata$Bradyrhizobium + 1))~log(mydata$Total_RL), data=mydata)
m4_1 = lm((log(mydata$Bradyrhizobium + 1))~log(mydata$TFRL2mm), data=mydata)
m4_2 = lm((log(mydata$Bradyrhizobium + 1))~log(mydata$FRL + 1), data=mydata)
m4_3 = lm((log(mydata$Bradyrhizobium + 1))~log(mydata$ExtraFRL), data=mydata)
m4_4 = lm((log(mydata$Bradyrhizobium + 1))~log(mydata$CoarseRL + 1), data=mydata)

#Summary
summary(m4)
summary(m4_1)
summary(m4_2)
summary(m4_3)
summary(m4_4)

## Flavobacterium

m5 = lm((log(mydata$Flavobacterium + 1))~log(mydata$Total_RL), data=mydata)
m5_1 = lm((log(mydata$Flavobacterium + 1))~log(mydata$TFRL2mm), data=mydata)
m5_2 = lm((log(mydata$Flavobacterium + 1))~log(mydata$FRL + 1), data=mydata)
m5_3 = lm((log(mydata$Flavobacterium + 1))~log(mydata$ExtraFRL), data=mydata)
m5_4 = lm((log(mydata$Flavobacterium + 1))~log(mydata$CoarseRL + 1), data=mydata)

#Summary
summary(m5)
summary(m5_1)
summary(m5_2)
summary(m5_3)
summary(m5_4)

## Gluconacetobacter

m6 = lm((log(mydata$Gluconacetobacter + 1))~log(mydata$Total_RL), data=mydata)
m6_1 = lm((log(mydata$Gluconacetobacter + 1))~log(mydata$TFRL2mm), data=mydata)
m6_2 = lm((log(mydata$Gluconacetobacter + 1))~log(mydata$FRL + 1), data=mydata)
m6_3 = lm((log(mydata$Gluconacetobacter + 1))~log(mydata$ExtraFRL), data=mydata)
m6_4 = lm((log(mydata$Gluconacetobacter + 1))~log(mydata$CoarseRL + 1), data=mydata)

#Summary
summary(m6)
summary(m6_1)
summary(m6_2)
summary(m6_3)
summary(m6_4)


## Lysobacter
m7 = lm((log(mydata$Lysobacter + 1))~log(mydata$Total_RL), data=mydata)
m7_1 = lm((log(mydata$Lysobacter + 1))~log(mydata$TFRL2mm), data=mydata)
m7_2 = lm((log(mydata$Lysobacter + 1))~log(mydata$FRL + 1), data=mydata)
m7_3 = lm((log(mydata$Lysobacter + 1))~log(mydata$ExtraFRL), data=mydata)
m7_4 = lm((log(mydata$Lysobacter + 1))~log(mydata$CoarseRL + 1), data=mydata)

#Summary
summary(m7)
summary(m7_1)
summary(m7_2)
summary(m7_3)
summary(m7_4)


## Microbacterium
m8 = lm((log(mydata$Microbacterium + 1))~log(mydata$Total_RL), data=mydata)
m8_1 = lm((log(mydata$Microbacterium + 1))~log(mydata$TFRL2mm), data=mydata)
m8_2 = lm((log(mydata$Microbacterium + 1))~log(mydata$FRL + 1), data=mydata)
m8_3 = lm((log(mydata$Microbacterium + 1))~log(mydata$ExtraFRL), data=mydata)
m8_4 = lm((log(mydata$Microbacterium + 1))~log(mydata$CoarseRL + 1), data=mydata)

#Summary
summary(m8)
summary(m8_1)
summary(m8_2)
summary(m8_3)
summary(m8_4)

## Paenibacillus
m9 = lm((log(mydata$Paenibacillus + 1))~log(mydata$Total_RL), data=mydata)
m9_1 = lm((log(mydata$Paenibacillus + 1))~log(mydata$TFRL2mm), data=mydata)
m9_2 = lm((log(mydata$Paenibacillus + 1))~log(mydata$FRL + 1), data=mydata)
m9_3 = lm((log(mydata$Paenibacillus + 1))~log(mydata$ExtraFRL), data=mydata)
m9_4 = lm((log(mydata$Paenibacillus + 1))~log(mydata$CoarseRL + 1), data=mydata)

#Summary
summary(m9)
summary(m9_1)
summary(m9_2)
summary(m9_3)
summary(m9_4)

## Pedobacter
m10 = lm((log(mydata$Pedobacter + 1))~log(mydata$Total_RL), data=mydata)
m10_1 = lm((log(mydata$Pedobacter + 1))~log(mydata$TFRL2mm), data=mydata)
m10_2 = lm((log(mydata$Pedobacter + 1))~log(mydata$FRL + 1), data=mydata)
m10_3 = lm((log(mydata$Pedobacter + 1))~log(mydata$ExtraFRL), data=mydata)
m10_4 = lm((log(mydata$Pedobacter + 1))~log(mydata$CoarseRL + 1), data=mydata)

#Summary
summary(m10)
summary(m10_1)
summary(m10_2)
summary(m10_3)
summary(m10_4)

## Phyllobacterium
m11 = lm((log(mydata$Phyllobacterium + 1))~log(mydata$Total_RL), data=mydata)
m11_1 = lm((log(mydata$Phyllobacterium + 1))~log(mydata$TFRL2mm), data=mydata)
m11_2 = lm((log(mydata$Phyllobacterium + 1))~log(mydata$FRL + 1), data=mydata)
m11_3 = lm((log(mydata$Phyllobacterium + 1))~log(mydata$ExtraFRL), data=mydata)
m11_4 = lm((log(mydata$Phyllobacterium + 1))~log(mydata$CoarseRL + 1), data=mydata)

#Summary
summary(m11)
summary(m11_1)
summary(m11_2)
summary(m11_3)
summary(m11_4)

## Pseudomonas
m12 = lm((log(mydata$Pseudomonas + 1))~log(mydata$Total_RL), data=mydata)
m12_1 = lm((log(mydata$Pseudomonas + 1))~log(mydata$TFRL2mm), data=mydata)
m12_2 = lm((log(mydata$Pseudomonas + 1))~log(mydata$FRL + 1), data=mydata)
m12_3 = lm((log(mydata$Pseudomonas + 1))~log(mydata$ExtraFRL), data=mydata)
m12_4 = lm((log(mydata$Pseudomonas + 1))~log(mydata$CoarseRL + 1), data=mydata)

#Summary
summary(m12)
summary(m12_1)
summary(m12_2)
summary(m12_3)
summary(m12_4)

## Rhizobium
m13 = lm((log(mydata$Rhizobium + 1))~log(mydata$Total_RL), data=mydata)
m13_1 = lm((log(mydata$Rhizobium + 1))~log(mydata$TFRL2mm), data=mydata)
m13_2 = lm((log(mydata$Rhizobium + 1))~log(mydata$FRL + 1), data=mydata)
m13_3 = lm((log(mydata$Rhizobium + 1))~log(mydata$ExtraFRL), data=mydata)
m13_4 = lm((log(mydata$Rhizobium + 1))~log(mydata$CoarseRL + 1), data=mydata)

#Summary
summary(m13)
summary(m13_1)
summary(m13_2)
summary(m13_3)
summary(m13_4)

## Rhodococcus
m14 = lm((log(mydata$Rhodococcus + 1))~log(mydata$Total_RL), data=mydata)
m14_1 = lm((log(mydata$Rhodococcus + 1))~log(mydata$TFRL2mm), data=mydata)
m14_2 = lm((log(mydata$Rhodococcus + 1))~log(mydata$FRL + 1), data=mydata)
m14_3 = lm((log(mydata$Rhodococcus + 1))~log(mydata$ExtraFRL), data=mydata)
m14_4 = lm((log(mydata$Rhodococcus + 1))~log(mydata$CoarseRL + 1), data=mydata)

#Summary
summary(m14)
summary(m14_1)
summary(m14_2)
summary(m14_3)
summary(m14_4)


## Rhodoplanes
m15 = lm((log(mydata$Rhodoplanes + 1))~log(mydata$Total_RL), data=mydata)
m15_1 = lm((log(mydata$Rhodoplanes + 1))~log(mydata$TFRL2mm), data=mydata)
m15_2 = lm((log(mydata$Rhodoplanes + 1))~log(mydata$FRL + 1), data=mydata)
m15_3 = lm((log(mydata$Rhodoplanes + 1))~log(mydata$ExtraFRL), data=mydata)
m15_4 = lm((log(mydata$Rhodoplanes + 1))~log(mydata$CoarseRL + 1), data=mydata)

#Summary
summary(m15)
summary(m15_1)
summary(m15_2)
summary(m15_3)
summary(m15_4)

## Skermanella
m16 = lm((log(mydata$Skermanella + 1))~log(mydata$Total_RL), data=mydata)
m16_1 = lm((log(mydata$Skermanella + 1))~log(mydata$TFRL2mm), data=mydata)
m16_2 = lm((log(mydata$Skermanella + 1))~log(mydata$FRL + 1), data=mydata)
m16_3 = lm((log(mydata$Skermanella + 1))~log(mydata$ExtraFRL), data=mydata)
m16_4 = lm((log(mydata$Skermanella + 1))~log(mydata$CoarseRL + 1), data=mydata)

#Summary
summary(m16)
summary(m16_1)
summary(m16_2)
summary(m16_3)
summary(m16_4)

##Stenotrophomonas
m17 = lm((log(mydata$Stenotrophomonas + 1))~log(mydata$Total_RL), data=mydata)
m17_1 = lm((log(mydata$Stenotrophomonas + 1))~log(mydata$TFRL2mm), data=mydata)
m17_2 = lm((log(mydata$Stenotrophomonas + 1))~log(mydata$FRL + 1), data=mydata)
m17_3 = lm((log(mydata$Stenotrophomonas + 1))~log(mydata$ExtraFRL), data=mydata)
m17_4 = lm((log(mydata$Stenotrophomonas + 1))~log(mydata$CoarseRL + 1), data=mydata)

#Summary
summary(m17)
summary(m17_1)
summary(m17_2)
summary(m17_3)
summary(m17_4)

## Variovorax
m18 = lm((log(mydata$Variovorax + 1))~log(mydata$Total_RL), data=mydata)
m18_1 = lm((log(mydata$Variovorax + 1))~log(mydata$TFRL2mm), data=mydata)
m18_2 = lm((log(mydata$Variovorax + 1))~log(mydata$FRL + 1), data=mydata)
m18_3 = lm((log(mydata$Variovorax + 1))~log(mydata$ExtraFRL), data=mydata)
m18_4 = lm((log(mydata$Variovorax + 1))~log(mydata$CoarseRL + 1), data=mydata)

#Summary
summary(m18)
summary(m18_1)
summary(m18_2)
summary(m18_3)
summary(m18_4)


## Xanthomonas
m19 = lm((log(mydata$Xanthomonas + 1))~log(mydata$Total_RL), data=mydata)
m19_1 = lm((log(mydata$Xanthomonas + 1))~log(mydata$TFRL2mm), data=mydata)
m19_2 = lm((log(mydata$Xanthomonas + 1))~log(mydata$FRL + 1), data=mydata)
m19_3 = lm((log(mydata$Xanthomonas + 1))~log(mydata$ExtraFRL), data=mydata)
m19_4 = lm((log(mydata$Xanthomonas + 1))~log(mydata$CoarseRL + 1), data=mydata)

#Summary
summary(m19)
summary(m19_1)
summary(m19_2)
summary(m19_3)
summary(m19_4)


#######################################

# Vegetative stage data
# Subset vegetative stage

vmydata = subset(mydata, Growth.Stage=="Vegetative")

#Fit linear models on log transformed bacterial genera abudnance and root traits: Vegetative

## Agrobacterium

m1 = lm((log(vmydata$Agrobacterium + 1))~log(vmydata$Total_RL), data=vmydata)
m1_1 = lm((log(vmydata$Agrobacterium + 1))~log(vmydata$TFRL2mm), data=vmydata)
m1_2 = lm((log(vmydata$Agrobacterium + 1))~log(vmydata$FRL + 1), data=vmydata)
m1_3 = lm((log(vmydata$Agrobacterium + 1))~log(vmydata$ExtraFRL), data=vmydata)
m1_4 = lm((log(vmydata$Agrobacterium + 1))~log(vmydata$CoarseRL + 1), data=vmydata)

#Summary
summary(m1)
summary(m1_1)
summary(m1_2)
summary(m1_3)
summary(m1_4)


## Arthrobacter

m2 = lm((log(vmydata$Arthrobacter + 1))~log(vmydata$Total_RL), data=vmydata)
m2_1 = lm((log(vmydata$Arthrobacter + 1))~log(vmydata$TFRL2mm), data=vmydata)
m2_2 = lm((log(vmydata$Arthrobacter + 1))~log(vmydata$FRL +1), data=vmydata)
m2_3 = lm((log(vmydata$Arthrobacter + 1))~log(vmydata$ExtraFRL), data=vmydata)
m2_4 = lm((log(vmydata$Arthrobacter + 1))~log(vmydata$CoarseRL + 1), data=vmydata)

#Summary
summary(m2)
summary(m2_1)
summary(m2_2)
summary(m2_3)
summary(m2_4)

## Bacillus
m3 = lm((log(vmydata$Bacillus + 1))~log(vmydata$Total_RL), data=vmydata)
m3_1 = lm((log(vmydata$Bacillus + 1))~log(vmydata$TFRL2mm), data=vmydata)
m3_2 = lm((log(vmydata$Bacillus + 1))~log(vmydata$FRL + 1), data=vmydata)
m3_3 = lm((log(vmydata$Bacillus + 1))~log(vmydata$ExtraFRL), data=vmydata)
m3_4 = lm((log(vmydata$Bacillus + 1))~log(vmydata$CoarseRL + 1), data=vmydata)

#Summary
summary(m3)
summary(m3_1)
summary(m3_2)
summary(m3_3)
summary(m3_4)

##Bradyrhizobium
m4 = lm((log(vmydata$Bradyrhizobium + 1))~log(vmydata$Total_RL), data=vmydata)
m4_1 = lm((log(vmydata$Bradyrhizobium + 1))~log(vmydata$TFRL2mm), data=vmydata)
m4_2 = lm((log(vmydata$Bradyrhizobium + 1))~log(vmydata$FRL + 1), data=vmydata)
m4_3 = lm((log(vmydata$Bradyrhizobium + 1))~log(vmydata$ExtraFRL), data=vmydata)
m4_4 = lm((log(vmydata$Bradyrhizobium + 1))~log(vmydata$CoarseRL + 1), data=vmydata)

#Summary
summary(m4)
summary(m4_1)
summary(m4_2)
summary(m4_3)
summary(m4_4)

## Flavobacterium

m5 = lm((log(vmydata$Flavobacterium + 1))~log(vmydata$Total_RL), data=vmydata)
m5_1 = lm((log(vmydata$Flavobacterium + 1))~log(vmydata$TFRL2mm), data=vmydata)
m5_2 = lm((log(vmydata$Flavobacterium + 1))~log(vmydata$FRL + 1), data=vmydata)
m5_3 = lm((log(vmydata$Flavobacterium + 1))~log(vmydata$ExtraFRL), data=vmydata)
m5_4 = lm((log(vmydata$Flavobacterium + 1))~log(vmydata$CoarseRL + 1), data=vmydata)

#Summary
summary(m5)
summary(m5_1)
summary(m5_2)
summary(m5_3)
summary(m5_4)

## Gluconacetobacter

m6 = lm((log(vmydata$Gluconacetobacter + 1))~log(vmydata$Total_RL), data=vmydata)
m6_1 = lm((log(vmydata$Gluconacetobacter + 1))~log(vmydata$TFRL2mm), data=vmydata)
m6_2 = lm((log(vmydata$Gluconacetobacter + 1))~log(vmydata$FRL + 1), data=vmydata)
m6_3 = lm((log(vmydata$Gluconacetobacter + 1))~log(vmydata$ExtraFRL), data=vmydata)
m6_4 = lm((log(vmydata$Gluconacetobacter + 1))~log(vmydata$CoarseRL + 1), data=vmydata)

#Summary
summary(m6)
summary(m6_1)
summary(m6_2)
summary(m6_3)
summary(m6_4)


## Lysobacter
m7 = lm((log(vmydata$Lysobacter + 1))~log(vmydata$Total_RL), data=vmydata)
m7_1 = lm((log(vmydata$Lysobacter + 1))~log(vmydata$TFRL2mm), data=vmydata)
m7_2 = lm((log(vmydata$Lysobacter + 1))~log(vmydata$FRL + 1), data=vmydata)
m7_3 = lm((log(vmydata$Lysobacter + 1))~log(vmydata$ExtraFRL), data=vmydata)
m7_4 = lm((log(vmydata$Lysobacter + 1))~log(vmydata$CoarseRL + 1), data=vmydata)

#Summary
summary(m7)
summary(m7_1)
summary(m7_2)
summary(m7_3)
summary(m7_4)


## Microbacterium
m8 = lm((log(vmydata$Microbacterium + 1))~log(vmydata$Total_RL), data=vmydata)
m8_1 = lm((log(vmydata$Microbacterium + 1))~log(vmydata$TFRL2mm), data=vmydata)
m8_2 = lm((log(vmydata$Microbacterium + 1))~log(vmydata$FRL + 1), data=vmydata)
m8_3 = lm((log(vmydata$Microbacterium + 1))~log(vmydata$ExtraFRL), data=vmydata)
m8_4 = lm((log(vmydata$Microbacterium + 1))~log(vmydata$CoarseRL + 1), data=vmydata)

#Summary
summary(m8)
summary(m8_1)
summary(m8_2)
summary(m8_3)
summary(m8_4)

## Paenibacillus
m9 = lm((log(vmydata$Paenibacillus + 1))~log(vmydata$Total_RL), data=vmydata)
m9_1 = lm((log(vmydata$Paenibacillus + 1))~log(vmydata$TFRL2mm), data=vmydata)
m9_2 = lm((log(vmydata$Paenibacillus + 1))~log(vmydata$FRL + 1), data=vmydata)
m9_3 = lm((log(vmydata$Paenibacillus + 1))~log(vmydata$ExtraFRL), data=vmydata)
m9_4 = lm((log(vmydata$Paenibacillus + 1))~log(vmydata$CoarseRL + 1), data=vmydata)

#Summary
summary(m9)
summary(m9_1)
summary(m9_2)
summary(m9_3)
summary(m9_4)

## Pedobacter
m10 = lm((log(vmydata$Pedobacter + 1))~log(vmydata$Total_RL), data=vmydata)
m10_1 = lm((log(vmydata$Pedobacter + 1))~log(vmydata$TFRL2mm), data=vmydata)
m10_2 = lm((log(vmydata$Pedobacter + 1))~log(vmydata$FRL + 1), data=vmydata)
m10_3 = lm((log(vmydata$Pedobacter + 1))~log(vmydata$ExtraFRL), data=vmydata)
m10_4 = lm((log(vmydata$Pedobacter + 1))~log(vmydata$CoarseRL + 1), data=vmydata)

#Summary
summary(m10)
summary(m10_1)
summary(m10_2)
summary(m10_3)
summary(m10_4)

## Phyllobacterium
m11 = lm((log(vmydata$Phyllobacterium + 1))~log(vmydata$Total_RL), data=vmydata)
m11_1 = lm((log(vmydata$Phyllobacterium + 1))~log(vmydata$TFRL2mm), data=vmydata)
m11_2 = lm((log(vmydata$Phyllobacterium + 1))~log(vmydata$FRL + 1), data=vmydata)
m11_3 = lm((log(vmydata$Phyllobacterium + 1))~log(vmydata$ExtraFRL), data=vmydata)
m11_4 = lm((log(vmydata$Phyllobacterium + 1))~log(vmydata$CoarseRL + 1), data=vmydata)

#Summary
summary(m11)
summary(m11_1)
summary(m11_2)
summary(m11_3)
summary(m11_4)

## Pseudomonas
m12 = lm((log(vmydata$Pseudomonas + 1))~log(vmydata$Total_RL), data=vmydata)
m12_1 = lm((log(vmydata$Pseudomonas + 1))~log(vmydata$TFRL2mm), data=vmydata)
m12_2 = lm((log(vmydata$Pseudomonas + 1))~log(vmydata$FRL + 1), data=vmydata)
m12_3 = lm((log(vmydata$Pseudomonas + 1))~log(vmydata$ExtraFRL), data=vmydata)
m12_4 = lm((log(vmydata$Pseudomonas + 1))~log(vmydata$CoarseRL + 1), data=vmydata)

#Summary
summary(m12)
summary(m12_1)
summary(m12_2)
summary(m12_3)
summary(m12_4)

## Rhizobium
m13 = lm((log(vmydata$Rhizobium + 1))~log(vmydata$Total_RL), data=vmydata)
m13_1 = lm((log(vmydata$Rhizobium + 1))~log(vmydata$TFRL2mm), data=vmydata)
m13_2 = lm((log(vmydata$Rhizobium + 1))~log(vmydata$FRL + 1), data=vmydata)
m13_3 = lm((log(vmydata$Rhizobium + 1))~log(vmydata$ExtraFRL), data=vmydata)
m13_4 = lm((log(vmydata$Rhizobium + 1))~log(vmydata$CoarseRL + 1), data=vmydata)

#Summary
summary(m13)
summary(m13_1)
summary(m13_2)
summary(m13_3)
summary(m13_4)

## Rhodococcus
m14 = lm((log(vmydata$Rhodococcus + 1))~log(vmydata$Total_RL), data=vmydata)
m14_1 = lm((log(vmydata$Rhodococcus + 1))~log(vmydata$TFRL2mm), data=vmydata)
m14_2 = lm((log(vmydata$Rhodococcus + 1))~log(vmydata$FRL + 1), data=vmydata)
m14_3 = lm((log(vmydata$Rhodococcus + 1))~log(vmydata$ExtraFRL), data=vmydata)
m14_4 = lm((log(vmydata$Rhodococcus + 1))~log(vmydata$CoarseRL + 1), data=vmydata)

#Summary
summary(m14)
summary(m14_1)
summary(m14_2)
summary(m14_3)
summary(m14_4)


## Rhodoplanes
m15 = lm((log(vmydata$Rhodoplanes + 1))~log(vmydata$Total_RL), data=vmydata)
m15_1 = lm((log(vmydata$Rhodoplanes + 1))~log(vmydata$TFRL2mm), data=vmydata)
m15_2 = lm((log(vmydata$Rhodoplanes + 1))~log(vmydata$FRL + 1), data=vmydata)
m15_3 = lm((log(vmydata$Rhodoplanes + 1))~log(vmydata$ExtraFRL), data=vmydata)
m15_4 = lm((log(vmydata$Rhodoplanes + 1))~log(vmydata$CoarseRL + 1), data=vmydata)

#Summary
summary(m15)
summary(m15_1)
summary(m15_2)
summary(m15_3)
summary(m15_4)

## Skermanella
m16 = lm((log(vmydata$Skermanella + 1))~log(vmydata$Total_RL), data=vmydata)
m16_1 = lm((log(vmydata$Skermanella + 1))~log(vmydata$TFRL2mm), data=vmydata)
m16_2 = lm((log(vmydata$Skermanella + 1))~log(vmydata$FRL + 1), data=vmydata)
m16_3 = lm((log(vmydata$Skermanella + 1))~log(vmydata$ExtraFRL), data=vmydata)
m16_4 = lm((log(vmydata$Skermanella + 1))~log(vmydata$CoarseRL), data=vmydata)

#Summary
summary(m16)
summary(m16_1)
summary(m16_2)
summary(m16_3)
summary(m16_4)

##Stenotrophomonas
m17 = lm((log(vmydata$Stenotrophomonas + 1))~log(vmydata$Total_RL), data=vmydata)
m17_1 = lm((log(vmydata$Stenotrophomonas + 1))~log(vmydata$TFRL2mm), data=vmydata)
m17_2 = lm((log(vmydata$Stenotrophomonas + 1))~log(vmydata$FRL + 1), data=vmydata)
m17_3 = lm((log(vmydata$Stenotrophomonas + 1))~log(vmydata$ExtraFRL), data=vmydata)
m17_4 = lm((log(vmydata$Stenotrophomonas + 1))~log(vmydata$CoarseRL + 1), data=vmydata)

#Summary
summary(m17)
summary(m17_1)
summary(m17_2)
summary(m17_3)
summary(m17_4)

## Variovorax
m18 = lm((log(vmydata$Variovorax + 1))~log(vmydata$Total_RL), data=vmydata)
m18_1 = lm((log(vmydata$Variovorax + 1))~log(vmydata$TFRL2mm), data=vmydata)
m18_2 = lm((log(vmydata$Variovorax + 1))~log(vmydata$FRL + 1), data=vmydata)
m18_3 = lm((log(vmydata$Variovorax + 1))~log(vmydata$ExtraFRL), data=vmydata)
m18_4 = lm((log(vmydata$Variovorax + 1))~log(vmydata$CoarseRL + 1), data=vmydata)

#Summary
summary(m18)
summary(m18_1)
summary(m18_2)
summary(m18_3)
summary(m18_4)


## Xanthomonas
m19 = lm((log(vmydata$Xanthomonas + 1))~log(vmydata$Total_RL), data=vmydata)
m19_1 = lm((log(vmydata$Xanthomonas + 1))~log(vmydata$TFRL2mm), data=vmydata)
m19_2 = lm((log(vmydata$Xanthomonas + 1))~log(vmydata$FRL + 1), data=vmydata)
m19_3 = lm((log(vmydata$Xanthomonas + 1))~log(vmydata$ExtraFRL), data=vmydata)
m19_4 = lm((log(vmydata$Xanthomonas + 1))~log(vmydata$CoarseRL + 1), data=vmydata)

#Summary
summary(m19)
summary(m19_1)
summary(m19_2)
summary(m19_3)
summary(m19_4)

## Additional Analysis: Revission #2: mixed effect models 

# Vegetative stage data
# Subset vegetative stage

# log transforming the bacterial abundances

vmydata = subset(mydata, Growth.Stage=="Vegetative")
vmydata$LogAgrobacterium = log(vmydata$Agrobacterium + 1)

vmydata$LogArthrobacter = log(vmydata$Arthrobacter + 1)


vmydata$LogBacillus =log(vmydata$Bacillus + 1)

vmydata$LogBradyrhizobium =log(vmydata$Bradyrhizobium + 1)

vmydata$LogFlavobacterium =log(vmydata$Flavobacterium + 1)

vmydata$LogGluconacetobacter =log(vmydata$Gluconacetobacter + 1)

vmydata$LogLysobacter =log(vmydata$Lysobacter + 1)

vmydata$LogMicrobacterium =log(vmydata$Microbacterium + 1)

vmydata$LogPaenibacillus =log(vmydata$Paenibacillus + 1)

vmydata$LogPedobacter =log(vmydata$Pedobacter + 1)

vmydata$LogPhyllobacterium =log(vmydata$Phyllobacterium + 1)

vmydata$LogPseudomonas =log(vmydata$Pseudomonas + 1)

vmydata$LogRhizobium =log(vmydata$Rhizobium + 1)

vmydata$LogRhodococcus =log(vmydata$Rhodococcus + 1)

vmydata$LogRhodoplanes =log(vmydata$Rhodoplanes + 1)

vmydata$LogSkermanella =log(vmydata$Skermanella + 1)

vmydata$LogStenotrophomonas =log(vmydata$Stenotrophomonas + 1)

vmydata$LogVariovorax =log(vmydata$Variovorax + 1)

vmydata$LogXanthomonas =log(vmydata$Xanthomonas + 1)

# Log transforming root traits


vmydata$LogTotal_RL = log(vmydata$Total_RL + 1)

vmydata$LogTFRL2mm = log(vmydata$TFRL2mm + 1)

vmydata$LogFRL = log(vmydata$FRL + 1)

vmydata$LogExtraFRL =log(vmydata$ExtraFRL + 1)

vmydata$LogCoarseRL = log(vmydata$CoarseRL + 1)

# make week as factor  or the mixed effect model

week = as.factor(vmydata$Week)
cline = as.factor(vmydata$Canola.Lines)

# Mixed efect models

## LogAgrobacterium

Agro_TRL_intc = lm(LogAgrobacterium ~ 1, data=vmydata)
Agro_TRL_min =  lm(LogAgrobacterium~ LogTotal_RL,data=vmydata)
Agro_TRL_full = lme(LogAgrobacterium~ LogTotal_RL , random=~1|week, method ="ML", data=vmydata)
Agro_TRL_canola = lme(LogAgrobacterium~ LogTotal_RL , random=list(~1|week,~1|cline), method ="ML", data=vmydata)
Agro_TRL_canola = lme(LogAgrobacterium~ LogTotal_RL , random=~1|cline, method ="ML", data=vmydata)


AIC(Agro_TRL_intc) #817.3006
AIC(Agro_TRL_min) # 815.6635
AIC(Agro_TRL_full) # 808.1027  
AIC(Agro_TRL_canola)# 809.1336


summary(Agro_TRL_intc)
summary(Agro_TRL_min)
summary(Agro_TRL_full)
summary(Agro_TRL_canola)

#
Agro_TRL_intc = lm(LogAgrobacterium ~ 1, data=vmydata)
Agro_TFRL_min =  lm(LogAgrobacterium~ LogTFRL2mm, data=vmydata)
Agro_TFRL_full = lme(LogAgrobacterium~ LogTFRL2mm , random=~1|week, method ="ML", data=vmydata)
Agro_TFRL_canola = lme(LogAgrobacterium~ LogTFRL2mm , random=~1|cline, method ="ML", data=vmydata)

AIC(Agro_TRL_intc) #817.3006
AIC(Agro_TFRL_min) # 816.2143
AIC(Agro_TFRL_full) # 807.9386  
AIC(Agro_TFRL_canola)# 818.2143

summary(Agro_TRL_intc)
summary(Agro_TFRL_min)
summary(Agro_TFRL_full)
summary(Agro_TFRL_canola)

#
Agro_TRL_intc = lm(LogAgrobacterium ~ 1, data=vmydata)
Agro_FRL_min =  lm(LogAgrobacterium~ LogFRL, data=vmydata)
Agro_FRL_full = lme(LogAgrobacterium~ LogFRL , random=~1|week, method ="ML", data=vmydata)
Agro_FRL_canola = lme(LogAgrobacterium~ LogFRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Agro_TRL_intc) #817.3006
AIC(Agro_FRL_min) # 810.1682
AIC(Agro_FRL_full) # 808.9524   
AIC(Agro_FRL_canola) #812.1424

summary(Agro_TRL_intc)
summary(Agro_FRL_min)
summary(Agro_FRL_full)
summary(Agro_FRL_canola)

#
Agro_TRL_intc = lm(LogAgrobacterium ~ 1, data=vmydata)
Agro_EFRL_min =  lm(LogAgrobacterium~ LogExtraFRL, data=vmydata)
Agro_EFRL_full = lme(LogAgrobacterium~ LogExtraFRL , random=~1|week, method ="ML", data=vmydata)
Agro_EFRL_canola = lme(LogAgrobacterium~ LogExtraFRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Agro_TRL_intc) #817.3006
AIC(Agro_EFRL_min) # 817.8206
AIC(Agro_EFRL_full) # 806.5879   
AIC(Agro_EFRL_canola)#819.8206

summary(Agro_TRL_intc)
summary(Agro_EFRL_min)
summary(Agro_EFRL_full)
summary(Agro_EFRL_canola)

#
Agro_TRL_intc = lm(LogAgrobacterium ~ 1, data=vmydata)
Agro_CRL_min =  lm(LogAgrobacterium~ LogCoarseRL, data=vmydata)
Agro_CRL_full = lme(LogAgrobacterium~ LogCoarseRL , random=~1|week, method ="ML", data=vmydata)
Agro_CRL_canola = lme(LogAgrobacterium~ LogCoarseRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Agro_TRL_intc) #817.3006
AIC(Agro_CRL_min) #  807.6826.
AIC(Agro_CRL_full)  # 808.6  
AIC(Agro_CRL_canola)# 809.6799

summary(Agro_TRL_intc)
summary(Agro_CRL_min)
summary(Agro_CRL_full)
summary(Agro_CRL_canola)

## LogArthrobacter
Arth_TRL_intc = lm(LogArthrobacter ~ 1, data=vmydata)
Arth_TRL_min =  lm(LogArthrobacter~ LogTotal_RL, data=vmydata)
Arth_TRL_full = lme(LogArthrobacter~ LogTotal_RL , random=~1|week, method ="ML", data=vmydata)
Arth_TRL_canola = lme(LogArthrobacter~ LogTotal_RL , random=~1|cline, method ="ML", data=vmydata)

AIC(Arth_TRL_intc) # 643.8118
AIC(Arth_TRL_min) # 632.2112
AIC(Arth_TRL_full) # 620.557
AIC(Arth_TRL_canola)#634.1522

summary(Arth_TRL_intc)
summary(Arth_TRL_min)
summary(Arth_TRL_full)
summary(Arth_TRL_canola)

#
Arth_TRL_intc = lm(LogArthrobacter ~ 1, data=vmydata)
Arth_TFRL_min =  lm(LogArthrobacter~ LogTFRL2mm, data=vmydata)
Arth_TFRL_full = lme(LogArthrobacter~ LogTFRL2mm , random=~1|week, method ="ML", data=vmydata)
Arth_TFRL_canola = lme(LogArthrobacter~ LogTFRL2mm , random=~1|cline, method ="ML", data=vmydata)

AIC(Arth_TRL_intc) # 643.8118
AIC(Arth_TFRL_min) # 633.5153
AIC(Arth_TFRL_full) # 620.5385   
AIC(Arth_TFRL_canola)# 635.4672

summary(Arth_TRL_intc)
summary(Arth_TFRL_min)
summary(Arth_TFRL_full)
summary(Arth_TFRL_canola)

#

Arth_TRL_intc = lm(LogArthrobacter ~ 1, data=vmydata)
Arth_FRL_min =  lm(LogArthrobacter~ LogFRL, data=vmydata)
Arth_FRL_full = lme(LogArthrobacter~ LogFRL , random=~1|week, method ="ML", data=vmydata)
Arth_FRL_canola = lme(LogArthrobacter~ LogFRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Arth_TRL_intc) # 643.8118
AIC(Arth_FRL_min) # 625.4933
AIC(Arth_FRL_full) # 619.3666   
AIC(Arth_FRL_canola)#627.256

summary(Arth_TRL_intc)
summary(Arth_FRL_min)
summary(Arth_FRL_full)
summary(Arth_FRL_canola)

#
Arth_TRL_intc = lm(LogArthrobacter ~ 1, data=vmydata)
Arth_EFRL_min =  lm(LogArthrobacter~ LogExtraFRL, data=vmydata)
Arth_EFRL_full = lme(LogArthrobacter~ LogExtraFRL , random=~1|week, method ="ML", data=vmydata)
Arth_EFRL_canola = lme(LogArthrobacter~ LogExtraFRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Arth_TRL_intc) # 643.8118
AIC(Arth_EFRL_min) # 636.9416
AIC(Arth_EFRL_full) # 620.2421   
AIC(Arth_EFRL_canola)# 638.9374

summary(Arth_TRL_intc)
summary(Arth_EFRL_min)
summary(Arth_EFRL_full)
summary(Arth_EFRL_canola)

#

Arth_TRL_intc = lm(LogArthrobacter ~ 1, data=vmydata)
Arth_CRL_min =  lm(LogArthrobacter~ LogCoarseRL, data=vmydata)
Arth_CRL_full = lme(LogArthrobacter~ LogCoarseRL , random=~1|week, method ="ML", data=vmydata)
Arth_CRL_canola = lme(LogArthrobacter~ LogCoarseRL , random=~1|cline, method ="ML", data=vmydata)


AIC(Arth_TRL_intc) # 643.8118
AIC(Arth_CRL_min) #  619.0473
AIC(Arth_CRL_full)  # 618.058 
AIC(Arth_CRL_canola)# 621.0473

summary(Arth_TRL_intc)
summary(Arth_CRL_min)
summary(Arth_CRL_full)
summary(Arth_CRL_canola)

# LogBacillus

Baci_TRL_intc = lm(LogBacillus ~ 1, data=vmydata)
Baci_TRL_min =  lm(LogBacillus~ LogTotal_RL, data=vmydata)
Baci_TRL_full = lme(LogBacillus~ LogTotal_RL , random=~1|week, method ="ML", data=vmydata)
Baci_TRL_canola = lme(LogBacillus~ LogTotal_RL , random=~1|cline, method ="ML", data=vmydata)

AIC(Baci_TRL_intc) # 765.0129
AIC(Baci_TRL_min) # 764.2373
AIC(Baci_TRL_full) # 766.185
AIC(Baci_TRL_canola) #766.2373

summary(Baci_TRL_intc)
summary(Baci_TRL_min)
summary(Baci_TRL_full)
summary(Baci_TRL_canola)

#
Baci_TRL_intc = lm(LogBacillus ~ 1, data=vmydata)
Baci_TFRL_min =  lm(LogBacillus~ LogTFRL2mm, data=vmydata)
Baci_TFRL_full = lme(LogBacillus~ LogTFRL2mm , random=~1|week, method ="ML", data=vmydata)
Baci_TFRL_canola = lme(LogBacillus~ LogTFRL2mm , random=~1|cline, method ="ML", data=vmydata)

AIC(Baci_TRL_intc) # 765.0129
AIC(Baci_TFRL_min) # 764.2985
AIC(Baci_TFRL_full) # 766.1877   
AIC(Baci_TFRL_canola) # 766.2985

summary(Baci_TRL_intc)
summary(Baci_TFRL_min)
summary(Baci_TFRL_full)
summary(Baci_TFRL_canola)

#

Baci_TRL_intc = lm(LogBacillus ~ 1, data=vmydata)
Baci_FRL_min =  lm(LogBacillus~ LogFRL, data=vmydata)
Baci_FRL_full = lme(LogBacillus~ LogFRL , random=~1|week, method ="ML", data=vmydata)
Baci_FRL_canola = lme(LogBacillus~ LogFRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Baci_TRL_intc) # 765.0129
AIC(Baci_FRL_min) # 763.1344
AIC(Baci_FRL_full) # 765.1344    
AIC(Baci_FRL_canola) # 765.1344

summary(Baci_TRL_intc)
summary(Baci_FRL_min)
summary(Baci_FRL_full)
summary(Baci_FRL_canola)

#
Baci_TRL_intc = lm(LogBacillus ~ 1, data=vmydata)
Baci_EFRL_min =  lm(LogBacillus~ LogExtraFRL, data=vmydata)
Baci_EFRL_full = lme(LogBacillus~ LogExtraFRL , random=~1|week, method ="ML", data=vmydata)
Baci_EFRL_canola = lme(LogBacillus~ LogExtraFRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Baci_TRL_intc) # 765.0129
AIC(Baci_EFRL_min) # 765.1629
AIC(Baci_EFRL_full) # 766.5274   
AIC(Baci_EFRL_canola) # 767.1629

summary(Baci_TRL_intc)
summary(Baci_EFRL_min)
summary(Baci_EFRL_full)
summary(Baci_EFRL_canola)

#

Baci_TRL_intc = lm(LogBacillus ~ 1, data=vmydata)
Baci_CRL_min =  lm(LogBacillus~ LogCoarseRL, data=vmydata)
Baci_CRL_full = lme(LogBacillus~ LogCoarseRL , random=~1|week, method ="ML", data=vmydata)
Baci_CRL_canola = lme(LogBacillus~ LogCoarseRL , random=~1|cline, method ="ML", data=vmydata)


AIC(Baci_TRL_intc) # 765.0129
AIC(Baci_CRL_min) #  763.3699
AIC(Baci_CRL_full)  # 765.3699 
AIC(Baci_CRL_canola) # 765.3699

summary(Baci_TRL_intc)
summary(Baci_CRL_min)
summary(Baci_CRL_full)
summary(Baci_CRL_canola)


# LogBradyrhizobium

Brady_TRL_intc = lm(LogBradyrhizobium ~ 1, data=vmydata)
Brady_TRL_min =  lm(LogBradyrhizobium~ LogTotal_RL, data=vmydata)
Brady_TRL_full = lme(LogBradyrhizobium~ LogTotal_RL , random=~1|week, method ="ML", data=vmydata)
Brady_TRL_canola = lme(LogBradyrhizobium~ LogTotal_RL , random=~1|cline, method ="ML", data=vmydata)

AIC(Brady_TRL_intc) # 636.9378
AIC(Brady_TRL_min) # 636.971
AIC(Brady_TRL_full) # 634.628
AIC(Brady_TRL_canola) # 638.971

summary(Brady_TRL_intc)
summary(Brady_TRL_min)
summary(Brady_TRL_full)
summary(Brady_TRL_canola)

#
Brady_TRL_intc = lm(LogBradyrhizobium ~ 1, data=vmydata)
Brady_TFRL_min =  lm(LogBradyrhizobium~ LogTFRL2mm, data=vmydata)
Brady_TFRL_full = lme(LogBradyrhizobium~ LogTFRL2mm , random=~1|week, method ="ML", data=vmydata)
Brady_TFRL_canola = lme(LogBradyrhizobium~ LogTFRL2mm , random=~1|cline, method ="ML", data=vmydata)

AIC(Brady_TRL_intc) # 636.9378
AIC(Brady_TFRL_min) # 637.2605
AIC(Brady_TFRL_full) # 634.5889  
AIC(Brady_TFRL_canola) # 639.2605

summary(Brady_TRL_intc)
summary(Brady_TFRL_min)
summary(Brady_TFRL_full)
summary(Brady_TFRL_canola)

#

Brady_TRL_intc = lm(LogBradyrhizobium ~ 1, data=vmydata)
Brady_FRL_min =  lm(LogBradyrhizobium~ LogFRL, data=vmydata)
Brady_FRL_full = lme(LogBradyrhizobium~ LogFRL , random=~1|week, method ="ML", data=vmydata)
Brady_FRL_canola = lme(LogBradyrhizobium~ LogFRL , random=~1|cline, method ="ML", data=vmydata)


AIC(Brady_TRL_intc) # 636.9378
AIC(Brady_FRL_min) # 635.842
AIC(Brady_FRL_full) # 635.0145   
AIC(Brady_FRL_canola) # 637.842

summary(Brady_TRL_intc)
summary(Brady_FRL_min)
summary(Brady_FRL_full)
summary(Brady_FRL_canola)

#
Brady_TRL_intc = lm(LogBradyrhizobium ~ 1, data=vmydata)
Brady_EFRL_min =  lm(LogBradyrhizobium~ LogExtraFRL, data=vmydata)
Brady_EFRL_full = lme(LogBradyrhizobium~ LogExtraFRL , random=~1|week, method ="ML", data=vmydata)
Brady_EFRL_canola = lme(LogBradyrhizobium~ LogExtraFRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Brady_TRL_intc) # 636.9378
AIC(Brady_EFRL_min) # 637.7426
AIC(Brady_EFRL_full) # 634.4478   
AIC(Brady_EFRL_canola)# 639.7426

summary(Brady_TRL_intc)
summary(Brady_EFRL_min)
summary(Brady_EFRL_full)
summary(Brady_EFRL_canola)

#

Brady_TRL_intc = lm(LogBradyrhizobium ~ 1, data=vmydata)
Brady_CRL_min =  lm(LogBradyrhizobium~ LogCoarseRL, data=vmydata)
Brady_CRL_full = lme(LogBradyrhizobium~ LogCoarseRL , random=~1|week, method ="ML", data=vmydata)
Brady_CRL_canola = lme(LogBradyrhizobium~ LogCoarseRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Brady_TRL_intc) # 636.9378
AIC(Brady_CRL_min) #  631.6471
AIC(Brady_CRL_full)  # 633.6471 
AIC(Brady_CRL_canola)# 633.6471

summary(Brady_TRL_intc)
summary(Brady_CRL_min)
summary(Brady_CRL_full)
summary(Brady_CRL_canola)


# LogFlavobacterium

Flavo_TRL_intc = lm(LogFlavobacterium ~ 1, data=vmydata)
Flavo_TRL_min =  lm(LogFlavobacterium~ LogTotal_RL, data=vmydata)
Flavo_TRL_full = lme(LogFlavobacterium~ LogTotal_RL , random=~1|week, method ="ML", data=vmydata)
Flavo_TRL_canola = lme(LogFlavobacterium~ LogTotal_RL , random=~1|cline, method ="ML", data=vmydata)

AIC(Flavo_TRL_intc) # 819.7526
AIC(Flavo_TRL_min) # 807.2047
AIC(Flavo_TRL_full) # 791.0682
AIC(Flavo_TRL_canola) # 809.2047

summary(Flavo_TRL_intc)
summary(Flavo_TRL_min)
summary(Flavo_TRL_full)
summary(Flavo_TRL_canola)

#
Flavo_TRL_intc = lm(LogFlavobacterium ~ 1, data=vmydata)
Flavo_TFRL_min =  lm(LogFlavobacterium~ LogTFRL2mm, data=vmydata)
Flavo_TFRL_full = lme(LogFlavobacterium~ LogTFRL2mm , random=~1|week, method ="ML", data=vmydata)
Flavo_TFRL_canola = lme(LogFlavobacterium~ LogTFRL2mm , random=~1|cline, method ="ML", data=vmydata)

AIC(Flavo_TRL_intc) # 819.7526
AIC(Flavo_TFRL_min) # 808.4615
AIC(Flavo_TFRL_full) # 791.0006   
AIC(Flavo_TFRL_canola)# 810.4615

summary(Flavo_TRL_intc)
summary(Flavo_TFRL_min)
summary(Flavo_TFRL_full)
summary(Flavo_TFRL_canola)

#

Flavo_TRL_intc = lm(LogFlavobacterium ~ 1, data=vmydata)
Flavo_FRL_min =  lm(LogFlavobacterium~ LogFRL, data=vmydata)
Flavo_FRL_full = lme(LogFlavobacterium~ LogFRL , random=~1|week, method ="ML", data=vmydata)
Flavo_FRL_canola = lme(LogFlavobacterium~ LogFRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Flavo_TRL_intc) # 819.7526
AIC(Flavo_FRL_min) # 802.107
AIC(Flavo_FRL_full) # 790.815    
AIC(Flavo_FRL_canola)# 804.107

summary(Flavo_TRL_intc)
summary(Flavo_FRL_min)
summary(Flavo_FRL_full)
summary(Flavo_FRL_canola)

#
Flavo_TRL_intc = lm(LogFlavobacterium ~ 1, data=vmydata)
Flavo_EFRL_min =  lm(LogFlavobacterium~ LogExtraFRL, data=vmydata)
Flavo_EFRL_full = lme(LogFlavobacterium~ LogExtraFRL , random=~1|week, method ="ML", data=vmydata)
Flavo_EFRL_canola = lme(LogFlavobacterium~ LogExtraFRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Flavo_TRL_intc) # 819.7526
AIC(Flavo_EFRL_min) # 811.3891
AIC(Flavo_EFRL_full) # 790.7954   
AIC(Flavo_EFRL_canola)# 813.3891

summary(Flavo_TRL_intc)
summary(Flavo_EFRL_min)
summary(Flavo_EFRL_full)
summary(Flavo_EFRL_canola)

#

Flavo_TRL_intc = lm(LogFlavobacterium ~ 1, data=vmydata)
Flavo_CRL_min =  lm(LogFlavobacterium~ LogCoarseRL, data=vmydata)
Flavo_CRL_full = lme(LogFlavobacterium~ LogCoarseRL , random=~1|week, method ="ML", data=vmydata)
Flavo_CRL_canola = lme(LogFlavobacterium~ LogCoarseRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Flavo_TRL_intc) # 819.7526
AIC(Flavo_CRL_min) #  802.0496
AIC(Flavo_CRL_full)  # 790.6967 
AIC(Flavo_CRL_canola) # 804.0496

summary(Flavo_TRL_intc)
summary(Flavo_CRL_min)
summary(Flavo_CRL_full)
summary(Flavo_CRL_canola)

# LogGluconacetobacter

Gluca_TRL_intc = lm(LogGluconacetobacter ~ 1, data=vmydata)
Gluca_TRL_min =  lm(LogGluconacetobacter~ LogTotal_RL, data=vmydata)
Gluca_TRL_full = lme(LogGluconacetobacter~ LogTotal_RL , random=~1|week, method ="ML", data=vmydata)
Gluca_TRL_canola = lme(LogGluconacetobacter~ LogTotal_RL , random=~1|cline, method ="ML", data=vmydata)

AIC(Gluca_TRL_intc) # 833.3612
AIC(Gluca_TRL_min) # 831.2576
AIC(Gluca_TRL_full) # 833.2576
AIC(Gluca_TRL_canola) # 833.1249

summary(Gluca_TRL_intc)
summary(Gluca_TRL_min)
summary(Gluca_TRL_full)
summary(Gluca_TRL_canola)

#
Gluca_TRL_intc = lm(LogGluconacetobacter ~ 1, data=vmydata)
Gluca_TFRL_min =  lm(LogGluconacetobacter~ LogTFRL2mm, data=vmydata)
Gluca_TFRL_full = lme(LogGluconacetobacter~ LogTFRL2mm , random=~1|week, method ="ML", data=vmydata)
Gluca_TFRL_canola = lme(LogGluconacetobacter~ LogTFRL2mm , random=~1|cline, method ="ML", data=vmydata)

AIC(Gluca_TRL_intc) # 833.3612
AIC(Gluca_TFRL_min) # 831.612
AIC(Gluca_TFRL_full) # 833.612   
AIC(Gluca_TFRL_canola)# 833.4815

summary(Gluca_TRL_intc)
summary(Gluca_TFRL_min)
summary(Gluca_TFRL_full)
summary(Gluca_TFRL_canola)

#

Gluca_TRL_intc = lm(LogGluconacetobacter ~ 1, data=vmydata)
Gluca_FRL_min =  lm(LogGluconacetobacter~ LogFRL, data=vmydata)
Gluca_FRL_full = lme(LogGluconacetobacter~ LogFRL , random=~1|week, method ="ML", data=vmydata)
Gluca_FRL_canola = lme(LogGluconacetobacter~ LogFRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Gluca_TRL_intc) # 833.3612
AIC(Gluca_FRL_min) # 828.8474
AIC(Gluca_FRL_full) # 830.8474    
AIC(Gluca_FRL_canola) # 830.6041

summary(Gluca_TRL_intc)
summary(Gluca_FRL_min)
summary(Gluca_FRL_full)
summary(Gluca_FRL_canola)

#
Gluca_TRL_intc = lm(LogGluconacetobacter ~ 1, data=vmydata)
Gluca_EFRL_min =  lm(LogGluconacetobacter~ LogExtraFRL, data=vmydata)
Gluca_EFRL_full = lme(LogGluconacetobacter~ LogExtraFRL , random=~1|week, method ="ML", data=vmydata)
Gluca_EFRL_canola = lme(LogGluconacetobacter~ LogExtraFRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Gluca_TRL_intc) # 833.3612
AIC(Gluca_EFRL_min) # 832.8365
AIC(Gluca_EFRL_full) # 834.8365   
AIC(Gluca_EFRL_canola) # 834.7453

summary(Gluca_TRL_intc)
summary(Gluca_EFRL_min)
summary(Gluca_EFRL_full)
summary(Gluca_EFRL_canola)

#

Gluca_TRL_intc = lm(LogGluconacetobacter ~ 1, data=vmydata)
Gluca_CRL_min =  lm(LogGluconacetobacter~ LogCoarseRL, data=vmydata)
Gluca_CRL_full = lme(LogGluconacetobacter~ LogCoarseRL , random=~1|week, method ="ML", data=vmydata)
Gluca_CRL_canola = lme(LogGluconacetobacter~ LogCoarseRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Gluca_TRL_intc) # 833.3612
AIC(Gluca_CRL_min) #  827.7682
AIC(Gluca_CRL_full)  # 829.7682 
AIC(Gluca_CRL_canola) # 829.6935

summary(Gluca_TRL_intc)
summary(Gluca_CRL_min)
summary(Gluca_CRL_full)
summary(Gluca_CRL_canola)


# LogLysobacter

Lysob_TRL_intc = lm(LogLysobacter ~ 1, data=vmydata)
Lysob_TRL_min =  lm(LogLysobacter~ LogTotal_RL, data=vmydata)
Lysob_TRL_full = lme(LogLysobacter~ LogTotal_RL , random=~1|week, method ="ML", data=vmydata)
Lysob_TRL_canola = lme(LogLysobacter~ LogTotal_RL , random=~1|cline, method ="ML", data=vmydata)

AIC(Lysob_TRL_intc) # 801.0342
AIC(Lysob_TRL_min) # 793.9852
AIC(Lysob_TRL_full) # 783.3093
AIC(Lysob_TRL_canola)# 795.9852

summary(Lysob_TRL_intc)
summary(Lysob_TRL_min)
summary(Lysob_TRL_full)
summary(Lysob_TRL_canola)

#
Lysob_TRL_intc = lm(LogLysobacter ~ 1, data=vmydata)
Lysob_TFRL_min =  lm(LogLysobacter~ LogTFRL2mm, data=vmydata)
Lysob_TFRL_full = lme(LogLysobacter~ LogTFRL2mm , random=~1|week, method ="ML", data=vmydata)
Lysob_TFRL_canola = lme(LogLysobacter~ LogTFRL2mm , random=~1|cline, method ="ML", data=vmydata)

AIC(Lysob_TRL_intc) # 801.0342
AIC(Lysob_TFRL_min) # 794.3988
AIC(Lysob_TFRL_full) # 783.3983   
AIC(Lysob_TFRL_canola) # 796.3988

summary(Lysob_TRL_intc)
summary(Lysob_TFRL_min)
summary(Lysob_TFRL_full)
summary(Lysob_TFRL_canola)

#

Lysob_TRL_intc = lm(LogLysobacter ~ 1, data=vmydata)
Lysob_FRL_min =  lm(LogLysobacter~ LogFRL, data=vmydata)
Lysob_FRL_full = lme(LogLysobacter~ LogFRL , random=~1|week, method ="ML", data=vmydata)
Lysob_FRL_canola = lme(LogLysobacter~ LogFRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Lysob_TRL_intc) # 801.0342
AIC(Lysob_FRL_min) # 790.2604
AIC(Lysob_FRL_full) # 782.384    
AIC(Lysob_FRL_canola) # 792.2604

summary(Lysob_TRL_intc)
summary(Lysob_FRL_min)
summary(Lysob_FRL_full)
summary(Lysob_FRL_canola)

#
Lysob_TRL_intc = lm(LogLysobacter ~ 1, data=vmydata)
Lysob_EFRL_min =  lm(LogLysobacter~ LogExtraFRL, data=vmydata)
Lysob_EFRL_full = lme(LogLysobacter~ LogExtraFRL , random=~1|week, method ="ML", data=vmydata)
Lysob_EFRL_canola = lme(LogLysobacter~ LogExtraFRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Lysob_TRL_intc) # 801.0342
AIC(Lysob_EFRL_min) # 796.6484
AIC(Lysob_EFRL_full) # 783.5287   
AIC(Lysob_EFRL_canola)# 798.6484

summary(Lysob_TRL_intc)
summary(Lysob_EFRL_min)
summary(Lysob_EFRL_full)
summary(Lysob_EFRL_canola)

#

Lysob_TRL_intc = lm(LogLysobacter ~ 1, data=vmydata)
Lysob_CRL_min =  lm(LogLysobacter~ LogCoarseRL, data=vmydata)
Lysob_CRL_full = lme(LogLysobacter~ LogCoarseRL , random=~1|week, method ="ML", data=vmydata)
Lysob_CRL_canola = lme(LogLysobacter~ LogCoarseRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Lysob_TRL_intc) # 801.0342
AIC(Lysob_CRL_min) #  792.4796
AIC(Lysob_CRL_full)  # 781.3166 
AIC(Lysob_CRL_canola) # 794.4796

summary(Lysob_TRL_intc)
summary(Lysob_CRL_min)
summary(Lysob_CRL_full)
summary(Lysob_CRL_canola)


# LogMicrobacterium


Microb_TRL_intc = lm(LogMicrobacterium ~ 1, data=vmydata)
Microb_TRL_min =  lm(LogMicrobacterium~ LogTotal_RL, data=vmydata)
Microb_TRL_full = lme(LogMicrobacterium~ LogTotal_RL , random=~1|week, method ="ML", data=vmydata)
Microb_TRL_canola = lme(LogMicrobacterium~ LogTotal_RL , random=~1|cline, method ="ML", data=vmydata)

AIC(Microb_TRL_intc) # 807.5441
AIC(Microb_TRL_min) # 805.9052
AIC(Microb_TRL_full) # 807.0883
AIC(Microb_TRL_canola)# 807.9052

summary(Microb_TRL_intc)
summary(Microb_TRL_min)
summary(Microb_TRL_full)
summary(Microb_TRL_canola)

#
Microb_TRL_intc = lm(LogMicrobacterium ~ 1, data=vmydata)
Microb_TFRL_min =  lm(LogMicrobacterium~ LogTFRL2mm, data=vmydata)
Microb_TFRL_full = lme(LogMicrobacterium~ LogTFRL2mm , random=~1|week, method ="ML", data=vmydata)
Microb_TFRL_canola = lme(LogMicrobacterium~ LogTFRL2mm , random=~1|cline, method ="ML", data=vmydata)

AIC(Microb_TRL_intc) # 807.5441
AIC(Microb_TFRL_min) # 806.151
AIC(Microb_TFRL_full) # 807.1182   
AIC(Microb_TFRL_canola) # 808.151

summary(Microb_TRL_intc)
summary(Microb_TFRL_min)
summary(Microb_TFRL_full)
summary(Microb_TFRL_canola)

#

Microb_TRL_intc = lm(LogMicrobacterium ~ 1, data=vmydata)
Microb_FRL_min =  lm(LogMicrobacterium~ LogFRL, data=vmydata)
Microb_FRL_full = lme(LogMicrobacterium~ LogFRL , random=~1|week, method ="ML", data=vmydata)
Microb_FRL_canola = lme(LogMicrobacterium~ LogFRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Microb_TRL_intc) # 807.5441
AIC(Microb_FRL_min) # 803.8496
AIC(Microb_FRL_full) # 805.6433    
AIC(Microb_FRL_canola) # 805.8496

summary(Microb_TRL_intc)
summary(Microb_FRL_min)
summary(Microb_FRL_full)
summary(Microb_FRL_canola)

#
Microb_TRL_intc = lm(LogMicrobacterium ~ 1, data=vmydata)
Microb_EFRL_min =  lm(LogMicrobacterium~ LogExtraFRL, data=vmydata)
Microb_EFRL_full = lme(LogMicrobacterium~ LogExtraFRL , random=~1|week, method ="ML", data=vmydata)
Microb_EFRL_canola = lme(LogMicrobacterium~ LogExtraFRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Microb_TRL_intc) # 807.5441
AIC(Microb_EFRL_min) # 806.9668
AIC(Microb_EFRL_full) # 807.3369   
AIC(Microb_EFRL_canola) # 808.9668

summary(Microb_TRL_intc)
summary(Microb_EFRL_min)
summary(Microb_EFRL_full)
summary(Microb_EFRL_canola)

#

Microb_TRL_intc = lm(LogMicrobacterium ~ 1, data=vmydata)
Microb_CRL_min =  lm(LogMicrobacterium~ LogCoarseRL, data=vmydata)
Microb_CRL_full = lme(LogMicrobacterium~ LogCoarseRL , random=~1|week, method ="ML", data=vmydata)
Microb_CRL_canola = lme(LogMicrobacterium~ LogCoarseRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Microb_TRL_intc) # 807.5441
AIC(Microb_CRL_min) #  804.3644
AIC(Microb_CRL_full)  # 806.3644 
AIC(Microb_CRL_canola) # 806.3644

summary(Microb_TRL_intc)
summary(Microb_CRL_min)
summary(Microb_CRL_full)
summary(Microb_CRL_canola)

# LogPaenibacillus

Paeni_TRL_intc = lm(LogPaenibacillus ~ 1, data=vmydata)
Paeni_TRL_min =  lm(LogPaenibacillus~ LogTotal_RL, data=vmydata)
Paeni_TRL_full = lme(LogPaenibacillus~ LogTotal_RL , random=~1|week, method ="ML", data=vmydata)
Paeni_TRL_canola = lme(LogPaenibacillus~ LogTotal_RL , random=~1|cline, method ="ML", data=vmydata)

AIC(Paeni_TRL_intc) # 714.7014
AIC(Paeni_TRL_min) # 711.0452
AIC(Paeni_TRL_full) # 713.0452
AIC(Paeni_TRL_canola) # 713.0452

summary(Paeni_TRL_intc)
summary(Paeni_TRL_min)
summary(Paeni_TRL_full)
summary(Paeni_TRL_canola)

#
Paeni_TRL_intc = lm(LogPaenibacillus ~ 1, data=vmydata)
Paeni_TFRL_min =  lm(LogPaenibacillus~ LogTFRL2mm, data=vmydata)
Paeni_TFRL_full = lme(LogPaenibacillus~ LogTFRL2mm , random=~1|week, method ="ML", data=vmydata)
Paeni_TFRL_canola = lme(LogPaenibacillus~ LogTFRL2mm , random=~1|cline, method ="ML", data=vmydata)

AIC(Paeni_TRL_intc) # 714.7014
AIC(Paeni_TFRL_min) # 711.3612
AIC(Paeni_TFRL_full) # 713.3612   
AIC(Paeni_TFRL_canola)# 713.3612

summary(Paeni_TRL_intc)
summary(Paeni_TFRL_min)
summary(Paeni_TFRL_full)
summary(Paeni_TFRL_canola)

#

Paeni_TRL_intc = lm(LogPaenibacillus ~ 1, data=vmydata)
Paeni_FRL_min =  lm(LogPaenibacillus~ LogFRL, data=vmydata)
Paeni_FRL_full = lme(LogPaenibacillus~ LogFRL , random=~1|week, method ="ML", data=vmydata)
Paeni_FRL_canola = lme(LogPaenibacillus~ LogFRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Paeni_TRL_intc) # 714.7014
AIC(Paeni_FRL_min) # 708.7719
AIC(Paeni_FRL_full) # 710.7719    
AIC(Paeni_FRL_canola) # 710.7719

summary(Paeni_TRL_intc)
summary(Paeni_FRL_min)
summary(Paeni_FRL_full)
summary(Paeni_FRL_canola)

#
Paeni_TRL_intc = lm(LogPaenibacillus ~ 1, data=vmydata)
Paeni_EFRL_min =  lm(LogPaenibacillus~ LogExtraFRL, data=vmydata)
Paeni_EFRL_full = lme(LogPaenibacillus~ LogExtraFRL , random=~1|week, method ="ML", data=vmydata)
Paeni_EFRL_canola = lme(LogPaenibacillus~ LogExtraFRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Paeni_TRL_intc) # 714.7014
AIC(Paeni_EFRL_min) # 712.9514
AIC(Paeni_EFRL_full) # 714.9514   
AIC(Paeni_EFRL_canola) # 714.9514

summary(Paeni_TRL_intc)
summary(Paeni_EFRL_min)
summary(Paeni_EFRL_full)
summary(Paeni_EFRL_canola)

#

Paeni_TRL_intc = lm(LogPaenibacillus ~ 1, data=vmydata)
Paeni_CRL_min =  lm(LogPaenibacillus~ LogCoarseRL, data=vmydata)
Paeni_CRL_full = lme(LogPaenibacillus~ LogCoarseRL , random=~1|week, method ="ML", data=vmydata)
Paeni_CRL_canola = lme(LogPaenibacillus~ LogCoarseRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Paeni_TRL_intc) # 714.7014
AIC(Paeni_CRL_min) #  709.2995
AIC(Paeni_CRL_full)  # 711.2995 
AIC(Paeni_CRL_canola) # 711.2995

summary(Paeni_TRL_intc)
summary(Paeni_CRL_min)
summary(Paeni_CRL_full)
summary(Paeni_CRL_canola)

# LogPedobacter

Pedob_TRL_intc = lm(LogPedobacter ~ 1, data=vmydata)
Pedob_TRL_min =  lm(LogPedobacter~ LogTotal_RL, data=vmydata)
Pedob_TRL_full = lme(LogPedobacter~ LogTotal_RL , random=~1|week, method ="ML", data=vmydata)
Pedob_TRL_canola = lme(LogPedobacter~ LogTotal_RL , random=~1|cline, method ="ML", data=vmydata)

AIC(Pedob_TRL_intc) # 834.6184
AIC(Pedob_TRL_min) # 791.4457
AIC(Pedob_TRL_full) # 773.7478
AIC(Pedob_TRL_canola) # 793.4457

summary(Pedob_TRL_intc)
summary(Pedob_TRL_min)
summary(Pedob_TRL_full)
summary(Pedob_TRL_canola)

#
Pedob_TRL_intc = lm(LogPedobacter ~ 1, data=vmydata)
Pedob_TFRL_min =  lm(LogPedobacter~ LogTFRL2mm, data=vmydata)
Pedob_TFRL_full = lme(LogPedobacter~ LogTFRL2mm , random=~1|week, method ="ML", data=vmydata)
Pedob_TFRL_canola = lme(LogPedobacter~ LogTFRL2mm , random=~1|cline, method ="ML", data=vmydata)

AIC(Pedob_TRL_intc) # 834.6184
AIC(Pedob_TFRL_min) # 792.838
AIC(Pedob_TFRL_full) # 773.8303   
AIC(Pedob_TFRL_canola) # 794.838

summary(Pedob_TRL_intc)
summary(Pedob_TFRL_min)
summary(Pedob_TFRL_full)
summary(Pedob_TFRL_canola)

#

Pedob_TRL_intc = lm(LogPedobacter ~ 1, data=vmydata)
Pedob_FRL_min =  lm(LogPedobacter~ LogFRL, data=vmydata)
Pedob_FRL_full = lme(LogPedobacter~ LogFRL , random=~1|week, method ="ML", data=vmydata)
Pedob_FRL_canola = lme(LogPedobacter~ LogFRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Pedob_TRL_intc) # 834.6184
AIC(Pedob_FRL_min) # 788.8233
AIC(Pedob_FRL_full) # 773.321    
AIC(Pedob_FRL_canola) # 790.8233

summary(Pedob_TRL_intc)
summary(Pedob_FRL_min)
summary(Pedob_FRL_full)
summary(Pedob_FRL_canola)

#
Pedob_TRL_intc = lm(LogPedobacter ~ 1, data=vmydata)
Pedob_EFRL_min =  lm(LogPedobacter~ LogExtraFRL, data=vmydata)
Pedob_EFRL_full = lme(LogPedobacter~ LogExtraFRL , random=~1|week, method ="ML", data=vmydata)
Pedob_EFRL_canola = lme(LogPedobacter~ LogExtraFRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Pedob_TRL_intc) # 834.6184
AIC(Pedob_EFRL_min) # 800.0128
AIC(Pedob_EFRL_full) # 775.4477   
AIC(Pedob_EFRL_canola) # 802.0128

summary(Pedob_TRL_intc)
summary(Pedob_EFRL_min)
summary(Pedob_EFRL_full)
summary(Pedob_EFRL_canola)

#

Pedob_TRL_intc = lm(LogPedobacter ~ 1, data=vmydata)
Pedob_CRL_min =  lm(LogPedobacter~ LogCoarseRL, data=vmydata)
Pedob_CRL_full = lme(LogPedobacter~ LogCoarseRL , random=~1|week, method ="ML", data=vmydata)
Pedob_CRL_canola = lme(LogPedobacter~ LogCoarseRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Pedob_TRL_intc) # 834.6184
AIC(Pedob_CRL_min) #  807.1873
AIC(Pedob_CRL_full)  # 779.3108 
AIC(Pedob_CRL_canola) # 809.1873

summary(Pedob_TRL_intc)
summary(Pedob_CRL_min)
summary(Pedob_CRL_full)
summary(Pedob_CRL_canola)

# LogPhyllobacterium

Phyll_TRL_intc = lm(LogPhyllobacterium ~ 1, data=vmydata)
Phyll_TRL_min =  lm(LogPhyllobacterium~ LogTotal_RL, data=vmydata)
Phyll_TRL_full = lme(LogPhyllobacterium~ LogTotal_RL , random=~1|week, method ="ML", data=vmydata)
Phyll_TRL_canola = lme(LogPhyllobacterium~ LogTotal_RL , random=~1|cline, method ="ML", data=vmydata)

AIC(Phyll_TRL_intc) # 792.9627
AIC(Phyll_TRL_min) # 791.4127
AIC(Phyll_TRL_full) # 793.3943
AIC(Phyll_TRL_canola)# 793.4127

summary(Phyll_TRL_intc)
summary(Phyll_TRL_min)
summary(Phyll_TRL_full)
summary(Phyll_TRL_canola)

#
Phyll_TRL_intc = lm(LogPhyllobacterium ~ 1, data=vmydata)
Phyll_TFRL_min =  lm(LogPhyllobacterium~ LogTFRL2mm, data=vmydata)
Phyll_TFRL_full = lme(LogPhyllobacterium~ LogTFRL2mm , random=~1|week, method ="ML", data=vmydata)
Phyll_TFRL_canola = lme(LogPhyllobacterium~ LogTFRL2mm , random=~1|cline, method ="ML", data=vmydata)

AIC(Phyll_TRL_intc) # 792.9627
AIC(Phyll_TFRL_min) # 791.0823
AIC(Phyll_TFRL_full) # 793.0769   
AIC(Phyll_TFRL_canola) # 793.0823

summary(Phyll_TRL_intc)
summary(Phyll_TFRL_min)
summary(Phyll_TFRL_full)
summary(Phyll_TFRL_canola)

#

Phyll_TRL_intc = lm(LogPhyllobacterium ~ 1, data=vmydata)
Phyll_FRL_min =  lm(LogPhyllobacterium~ LogFRL, data=vmydata)
Phyll_FRL_full = lme(LogPhyllobacterium~ LogFRL , random=~1|week, method ="ML", data=vmydata)
Phyll_FRL_canola = lme(LogPhyllobacterium~ LogFRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Phyll_TRL_intc) # 792.9627
AIC(Phyll_FRL_min) # 793.6982
AIC(Phyll_FRL_full) # 795.687    
AIC(Phyll_FRL_canola) # 795.6982

summary(Phyll_TRL_intc)
summary(Phyll_FRL_min)
summary(Phyll_FRL_full)
summary(Phyll_FRL_canola)

#
Phyll_TRL_intc = lm(LogPhyllobacterium ~ 1, data=vmydata)
Phyll_EFRL_min =  lm(LogPhyllobacterium~ LogExtraFRL, data=vmydata)
Phyll_EFRL_full = lme(LogPhyllobacterium~ LogExtraFRL , random=~1|week, method ="ML", data=vmydata)
Phyll_EFRL_canola = lme(LogPhyllobacterium~ LogExtraFRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Phyll_TRL_intc) # 792.9627
AIC(Phyll_EFRL_min) # 790.2013
AIC(Phyll_EFRL_full) # 775.4477   
AIC(Phyll_EFRL_canola)# 792.2013

summary(Phyll_TRL_intc)
summary(Phyll_EFRL_min)
summary(Phyll_EFRL_full)
summary(Phyll_EFRL_canola)

#

Phyll_TRL_intc = lm(LogPhyllobacterium ~ 1, data=vmydata)
Phyll_CRL_min =  lm(LogPhyllobacterium~ LogCoarseRL, data=vmydata)
Phyll_CRL_full = lme(LogPhyllobacterium~ LogCoarseRL , random=~1|week, method ="ML", data=vmydata)
Phyll_CRL_canola = lme(LogPhyllobacterium~ LogCoarseRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Phyll_TRL_intc) # 792.9627
AIC(Phyll_CRL_min) #  794.8905
AIC(Phyll_CRL_full)  # 779.3108 
AIC(Phyll_CRL_canola) # 796.8905

summary(Phyll_TRL_intc)
summary(Phyll_CRL_min)
summary(Phyll_CRL_full)
summary(Phyll_CRL_canola)

# LogPseudomonas


Pseud_TRL_intc = lm(LogPseudomonas ~ 1, data=vmydata)
Pseud_TRL_min =  lm(LogPseudomonas~ LogTotal_RL, data=vmydata)
Pseud_TRL_full = lme(LogPseudomonas~ LogTotal_RL , random=~1|week, method ="ML", data=vmydata)
Pseud_TRL_canola = lme(LogPseudomonas~ LogTotal_RL , random=~1|cline, method ="ML", data=vmydata)

AIC(Pseud_TRL_intc) # 828.8763
AIC(Pseud_TRL_min) # 814.1951
AIC(Pseud_TRL_full) # 804.3425
AIC(Pseud_TRL_canola) # 816.1951

summary(Pseud_TRL_intc)
summary(Pseud_TRL_min)
summary(Pseud_TRL_full)
summary(Pseud_TRL_canola)

#
Pseud_TRL_intc = lm(LogPseudomonas ~ 1, data=vmydata)
Pseud_TFRL_min =  lm(LogPseudomonas~ LogTFRL2mm, data=vmydata)
Pseud_TFRL_full = lme(LogPseudomonas~ LogTFRL2mm , random=~1|week, method ="ML", data=vmydata)
Pseud_TFRL_canola = lme(LogPseudomonas~ LogTFRL2mm , random=~1|cline, method ="ML", data=vmydata)

AIC(Pseud_TRL_intc) # 828.8763
AIC(Pseud_TFRL_min) # 815.6096
AIC(Pseud_TFRL_full) # 804.4157   
AIC(Pseud_TFRL_canola) # 817.6096

summary(Pseud_TRL_intc)
summary(Pseud_TFRL_min)
summary(Pseud_TFRL_full)
summary(Pseud_TFRL_canola)

#

Pseud_TRL_intc = lm(LogPseudomonas ~ 1, data=vmydata)
Pseud_FRL_min =  lm(LogPseudomonas~ LogFRL, data=vmydata)
Pseud_FRL_full = lme(LogPseudomonas~ LogFRL , random=~1|week, method ="ML", data=vmydata)
Pseud_FRL_canola = lme(LogPseudomonas~ LogFRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Pseud_TRL_intc) # 828.8763
AIC(Pseud_FRL_min) # 814.1481
AIC(Pseud_FRL_full) # 804.4037    
AIC(Pseud_FRL_canola) # 816.1481

summary(Pseud_TRL_intc)
summary(Pseud_FRL_min)
summary(Pseud_FRL_full)
summary(Pseud_FRL_canola)

#
Pseud_TRL_intc = lm(LogPseudomonas ~ 1, data=vmydata)
Pseud_EFRL_min =  lm(LogPseudomonas~ LogExtraFRL, data=vmydata)
Pseud_EFRL_full = lme(LogPseudomonas~ LogExtraFRL , random=~1|week, method ="ML", data=vmydata)
Pseud_EFRL_canola = lme(LogPseudomonas~ LogExtraFRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Pseud_TRL_intc) # 828.8763
AIC(Pseud_EFRL_min) # 817.4361
AIC(Pseud_EFRL_full) # 804.4231   
AIC(Pseud_EFRL_canola) # 819.4361

summary(Pseud_TRL_intc)
summary(Pseud_EFRL_min)
summary(Pseud_EFRL_full)
summary(Pseud_EFRL_canola)

#

Pseud_TRL_intc = lm(LogPseudomonas ~ 1, data=vmydata)
Pseud_CRL_min =  lm(LogPseudomonas~ LogCoarseRL, data=vmydata)
Pseud_CRL_full = lme(LogPseudomonas~ LogCoarseRL , random=~1|week, method ="ML", data=vmydata)
Pseud_CRL_canola = lme(LogPseudomonas~ LogCoarseRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Pseud_TRL_intc) # 828.8763
AIC(Pseud_CRL_min) #  807.5735
AIC(Pseud_CRL_full)  # 803.2946 
AIC(Pseud_CRL_canola)# 809.5735

summary(Pseud_TRL_intc)
summary(Pseud_CRL_min)
summary(Pseud_CRL_full)
summary(Pseud_CRL_canola)


# LogRhizobium

Rhizo_TRL_intc = lm(LogRhizobium ~ 1, data=vmydata)
Rhizo_TRL_min =  lm(LogRhizobium~ LogTotal_RL, data=vmydata)
Rhizo_TRL_full = lme(LogRhizobium~ LogTotal_RL , random=~1|week, method ="ML", data=vmydata)
Rhizo_TRL_canola = lme(LogRhizobium~ LogTotal_RL , random=~1|cline, method ="ML", data=vmydata)

AIC(Rhizo_TRL_intc) # 815.1382
AIC(Rhizo_TRL_min) # 809.3474
AIC(Rhizo_TRL_full) # 809.9591
AIC(Rhizo_TRL_canola) # 810.2205

summary(Rhizo_TRL_intc)
summary(Rhizo_TRL_min)
summary(Rhizo_TRL_full)
summary(Rhizo_TRL_canola)

#
Rhizo_TRL_intc = lm(LogRhizobium ~ 1, data=vmydata)
Rhizo_TFRL_min =  lm(LogRhizobium~ LogTFRL2mm, data=vmydata)
Rhizo_TFRL_full = lme(LogRhizobium~ LogTFRL2mm , random=~1|week, method ="ML", data=vmydata)
Rhizo_TFRL_canola = lme(LogRhizobium~ LogTFRL2mm , random=~1|cline, method ="ML", data=vmydata)

AIC(Rhizo_TRL_intc) # 815.1382
AIC(Rhizo_TFRL_min) # 809.8495
AIC(Rhizo_TFRL_full) # 810.1168   
AIC(Rhizo_TFRL_canola) # 810.7412

summary(Rhizo_TRL_intc)
summary(Rhizo_TFRL_min)
summary(Rhizo_TFRL_full)
summary(Rhizo_TFRL_canola)

#

Rhizo_TRL_intc = lm(LogRhizobium ~ 1, data=vmydata)
Rhizo_FRL_min =  lm(LogRhizobium~ LogFRL, data=vmydata)
Rhizo_FRL_full = lme(LogRhizobium~ LogFRL , random=~1|week, method ="ML", data=vmydata)
Rhizo_FRL_canola = lme(LogRhizobium~ LogFRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Rhizo_TRL_intc) # 815.1382
AIC(Rhizo_FRL_min) # 806.6431
AIC(Rhizo_FRL_full) # 808.4139    
AIC(Rhizo_FRL_canola) # 807.3254

summary(Rhizo_TRL_intc)
summary(Rhizo_FRL_min)
summary(Rhizo_FRL_full)
summary(Rhizo_FRL_canola)

#
Rhizo_TRL_intc = lm(LogRhizobium ~ 1, data=vmydata)
Rhizo_EFRL_min =  lm(LogRhizobium~ LogExtraFRL, data=vmydata)
Rhizo_EFRL_full = lme(LogRhizobium~ LogExtraFRL , random=~1|week, method ="ML", data=vmydata)
Rhizo_EFRL_canola = lme(LogRhizobium~ LogExtraFRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Rhizo_TRL_intc) # 815.1382
AIC(Rhizo_EFRL_min) # 811.2033
AIC(Rhizo_EFRL_full) # 810.418   
AIC(Rhizo_EFRL_canola) # 812.279

summary(Rhizo_TRL_intc)
summary(Rhizo_EFRL_min)
summary(Rhizo_EFRL_full)
summary(Rhizo_EFRL_canola)

#

Rhizo_TRL_intc = lm(LogRhizobium ~ 1, data=vmydata)
Rhizo_CRL_min =  lm(LogRhizobium~ LogCoarseRL, data=vmydata)
Rhizo_CRL_full = lme(LogRhizobium~ LogCoarseRL , random=~1|week, method ="ML", data=vmydata)
Rhizo_CRL_canola = lme(LogRhizobium~ LogCoarseRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Rhizo_TRL_intc) # 815.1382
AIC(Rhizo_CRL_min) #  806.2488
AIC(Rhizo_CRL_full)  # 807.7289 
AIC(Rhizo_CRL_canola) # 807.4049

summary(Rhizo_TRL_intc)
summary(Rhizo_CRL_min)
summary(Rhizo_CRL_full)
summary(Rhizo_CRL_canola)

# LogRhodococcus

Rhodo_TRL_intc = lm(LogRhodococcus ~ 1, data=vmydata)
Rhodo_TRL_min =  lm(LogRhodococcus~ LogTotal_RL, data=vmydata)
Rhodo_TRL_full = lme(LogRhodococcus~ LogTotal_RL , random=~1|week, method ="ML", data=vmydata)
Rhodo_TRL_canola = lme(LogRhodococcus~ LogTotal_RL , random=~1|cline, method ="ML", data=vmydata)

AIC(Rhodo_TRL_intc) # 695.5589
AIC(Rhodo_TRL_min) # 697.4227
AIC(Rhodo_TRL_full) # 696.1386
AIC(Rhodo_TRL_canola)# 699.4227

summary(Rhodo_TRL_intc)
summary(Rhodo_TRL_min)
summary(Rhodo_TRL_full)
summary(Rhodo_TRL_canola)

#
Rhodo_TRL_intc = lm(LogRhodococcus ~ 1, data=vmydata)
Rhodo_TFRL_min =  lm(LogRhodococcus~ LogTFRL2mm, data=vmydata)
Rhodo_TFRL_full = lme(LogRhodococcus~ LogTFRL2mm , random=~1|week, method ="ML", data=vmydata)
Rhodo_TFRL_canola = lme(LogRhodococcus~ LogTFRL2mm , random=~1|cline, method ="ML", data=vmydata)

AIC(Rhodo_TRL_intc) # 695.5589
AIC(Rhodo_TFRL_min) # 697.446
AIC(Rhodo_TFRL_full) # 696.2856   
AIC(Rhodo_TFRL_canola)# 699.446

summary(Rhodo_TRL_intc)
summary(Rhodo_TFRL_min)
summary(Rhodo_TFRL_full)
summary(Rhodo_TFRL_canola)

#

Rhodo_TRL_intc = lm(LogRhodococcus ~ 1, data=vmydata)
Rhodo_FRL_min =  lm(LogRhodococcus~ LogFRL, data=vmydata)
Rhodo_FRL_full = lme(LogRhodococcus~ LogFRL , random=~1|week, method ="ML", data=vmydata)
Rhodo_FRL_canola = lme(LogRhodococcus~ LogFRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Rhodo_TRL_intc) # 695.5589
AIC(Rhodo_FRL_min) # 697.0722
AIC(Rhodo_FRL_full) # 697.0549    
AIC(Rhodo_FRL_canola) # 699.0722

summary(Rhodo_TRL_intc)
summary(Rhodo_FRL_min)
summary(Rhodo_FRL_full)
summary(Rhodo_FRL_canola)

#
Rhodo_TRL_intc = lm(LogRhodococcus ~ 1, data=vmydata)
Rhodo_EFRL_min =  lm(LogRhodococcus~ LogExtraFRL, data=vmydata)
Rhodo_EFRL_full = lme(LogRhodococcus~ LogExtraFRL , random=~1|week, method ="ML", data=vmydata)
Rhodo_EFRL_canola = lme(LogRhodococcus~ LogExtraFRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Rhodo_TRL_intc) # 695.5589
AIC(Rhodo_EFRL_min) # 697.5141
AIC(Rhodo_EFRL_full) # 696.2347   
AIC(Rhodo_EFRL_canola) # 699.5141

summary(Rhodo_TRL_intc)
summary(Rhodo_EFRL_min)
summary(Rhodo_EFRL_full)
summary(Rhodo_EFRL_canola)

#

Rhodo_TRL_intc = lm(LogRhodococcus ~ 1, data=vmydata)
Rhodo_CRL_min =  lm(LogRhodococcus~ LogCoarseRL, data=vmydata)
Rhodo_CRL_full = lme(LogRhodococcus~ LogCoarseRL , random=~1|week, method ="ML", data=vmydata)
Rhodo_CRL_canola = lme(LogRhodococcus~ LogCoarseRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Rhodo_TRL_intc) # 695.5589
AIC(Rhodo_CRL_min) #  696.8365
AIC(Rhodo_CRL_full)  # 696.3075 
AIC(Rhodo_CRL_canola)# 698.8365

summary(Rhodo_TRL_intc)
summary(Rhodo_CRL_min)
summary(Rhodo_CRL_full)
summary(Rhodo_CRL_canola)

# LogRhodoplanes


Rhodop_TRL_intc = lm(LogRhodoplanes ~ 1, data=vmydata)
Rhodop_TRL_min =  lm(LogRhodoplanes~ LogTotal_RL, data=vmydata)
Rhodop_TRL_full = lme(LogRhodoplanes~ LogTotal_RL , random=~1|week, method ="ML", data=vmydata)
Rhodop_TRL_canola = lme(LogRhodoplanes~ LogTotal_RL , random=~1|cline, method ="ML", data=vmydata)

AIC(Rhodop_TRL_intc) # 833.0205
AIC(Rhodop_TRL_min) # 833.2084
AIC(Rhodop_TRL_full) # 835.2084
AIC(Rhodop_TRL_canola)# 835.2084

summary(Rhodop_TRL_intc)
summary(Rhodop_TRL_min)
summary(Rhodop_TRL_full)
summary(Rhodop_TRL_canola)

#
Rhodop_TRL_intc = lm(LogRhodoplanes ~ 1, data=vmydata)
Rhodop_TFRL_min =  lm(LogRhodoplanes~ LogTFRL2mm, data=vmydata)
Rhodop_TFRL_full = lme(LogRhodoplanes~ LogTFRL2mm , random=~1|week, method ="ML", data=vmydata)
Rhodop_TFRL_canola = lme(LogRhodoplanes~ LogTFRL2mm , random=~1|cline, method ="ML", data=vmydata)

AIC(Rhodop_TRL_intc) # 833.0205
AIC(Rhodop_TFRL_min) # 833.2867
AIC(Rhodop_TFRL_full) # 835.2867   
AIC(Rhodop_TFRL_canola)# 835.2867

summary(Rhodop_TRL_intc)
summary(Rhodop_TFRL_min)
summary(Rhodop_TFRL_full)
summary(Rhodop_TFRL_canola)

#

Rhodop_TRL_intc = lm(LogRhodoplanes ~ 1, data=vmydata)
Rhodop_FRL_min =  lm(LogRhodoplanes~ LogFRL, data=vmydata)
Rhodop_FRL_full = lme(LogRhodoplanes~ LogFRL , random=~1|week, method ="ML", data=vmydata)
Rhodop_FRL_canola = lme(LogRhodoplanes~ LogFRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Rhodop_TRL_intc) # 833.0205
AIC(Rhodop_FRL_min) # 834.7588
AIC(Rhodop_FRL_full) # 836.7588    
AIC(Rhodop_FRL_canola) # 836.7588

summary(Rhodop_TRL_intc)
summary(Rhodop_FRL_min)
summary(Rhodop_FRL_full)
summary(Rhodop_FRL_canola)

#
Rhodop_TRL_intc = lm(LogRhodoplanes ~ 1, data=vmydata)
Rhodop_EFRL_min =  lm(LogRhodoplanes~ LogExtraFRL, data=vmydata)
Rhodop_EFRL_full = lme(LogRhodoplanes~ LogExtraFRL , random=~1|week, method ="ML", data=vmydata)
Rhodop_EFRL_canola = lme(LogRhodoplanes~ LogExtraFRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Rhodop_TRL_intc) # 833.0205
AIC(Rhodop_EFRL_min) # 832.641
AIC(Rhodop_EFRL_full) # 834.641   
AIC(Rhodop_EFRL_canola)# 834.641

summary(Rhodop_TRL_intc)
summary(Rhodop_EFRL_min)
summary(Rhodop_EFRL_full)
summary(Rhodop_EFRL_canola)

#

Rhodop_TRL_intc = lm(LogRhodoplanes ~ 1, data=vmydata)
Rhodop_CRL_min =  lm(LogRhodoplanes~ LogCoarseRL, data=vmydata)
Rhodop_CRL_full = lme(LogRhodoplanes~ LogCoarseRL , random=~1|week, method ="ML", data=vmydata)
Rhodop_CRL_canola = lme(LogRhodoplanes~ LogCoarseRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Rhodop_TRL_intc) # 833.0205
AIC(Rhodop_CRL_min) #  833.6836
AIC(Rhodop_CRL_full)  # 835.6836 
AIC(Rhodop_CRL_canola) # 835.6836

summary(Rhodop_TRL_intc)
summary(Rhodop_CRL_min)
summary(Rhodop_CRL_full)
summary(Rhodop_CRL_canola)


# LogSkermanella


Skerma_TRL_intc = lm(LogSkermanella ~ 1, data=vmydata)
Skerma_TRL_min =  lm(LogSkermanella~ LogTotal_RL, data=vmydata)
Skerma_TRL_full = lme(LogSkermanella~ LogTotal_RL , random=~1|week, method ="ML", data=vmydata)
Skerma_TRL_canola = lme(LogSkermanella~ LogTotal_RL , random=~1|cline, method ="ML", data=vmydata)

AIC(Skerma_TRL_intc) # 691.7127
AIC(Skerma_TRL_min) # 691.8137
AIC(Skerma_TRL_full) # 691.6617
AIC(Skerma_TRL_canola)# 693.8137

summary(Skerma_TRL_intc)
summary(Skerma_TRL_min)
summary(Skerma_TRL_full)
summary(Skerma_TRL_canola)

#
Skerma_TRL_intc = lm(LogSkermanella ~ 1, data=vmydata)
Skerma_TFRL_min =  lm(LogSkermanella~ LogTFRL2mm, data=vmydata)
Skerma_TFRL_full = lme(LogSkermanella~ LogTFRL2mm , random=~1|week, method ="ML", data=vmydata)
Skerma_TFRL_canola = lme(LogSkermanella~ LogTFRL2mm , random=~1|cline, method ="ML", data=vmydata)


AIC(Skerma_TRL_intc) # 691.7127
AIC(Skerma_TFRL_min) # 692.0873
AIC(Skerma_TFRL_full) # 691.6073   
AIC(Skerma_TFRL_canola)# 694.0873

summary(Skerma_TRL_intc)
summary(Skerma_TFRL_min)
summary(Skerma_TFRL_full)
summary(Skerma_TFRL_canola)

#

Skerma_TRL_intc = lm(LogSkermanella ~ 1, data=vmydata)
Skerma_FRL_min =  lm(LogSkermanella~ LogFRL, data=vmydata)
Skerma_FRL_full = lme(LogSkermanella~ LogFRL , random=~1|week, method ="ML", data=vmydata)
Skerma_FRL_canola = lme(LogSkermanella~ LogFRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Skerma_TRL_intc) # 691.7127
AIC(Skerma_FRL_min) # 686.0122
AIC(Skerma_FRL_full) # 688.0122    
AIC(Skerma_FRL_canola)# 688.0122

summary(Skerma_TRL_intc)
summary(Skerma_FRL_min)
summary(Skerma_FRL_full)
summary(Skerma_FRL_canola)

#
Skerma_TRL_intc = lm(LogSkermanella ~ 1, data=vmydata)
Skerma_EFRL_min =  lm(LogSkermanella~ LogExtraFRL, data=vmydata)
Skerma_EFRL_full = lme(LogSkermanella~ LogExtraFRL , random=~1|week, method ="ML", data=vmydata)
Skerma_EFRL_canola = lme(LogSkermanella~ LogExtraFRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Skerma_TRL_intc) # 691.7127
AIC(Skerma_EFRL_min) # 693.1454
AIC(Skerma_EFRL_full) # 690.946   
AIC(Skerma_EFRL_canola) # 695.1454

summary(Skerma_TRL_intc)
summary(Skerma_EFRL_min)
summary(Skerma_EFRL_full)
summary(Skerma_EFRL_canola)

#

Skerma_TRL_intc = lm(LogSkermanella ~ 1, data=vmydata)
Skerma_CRL_min =  lm(LogSkermanella~ LogCoarseRL, data=vmydata)
Skerma_CRL_full = lme(LogSkermanella~ LogCoarseRL , random=~1|week, method ="ML", data=vmydata)
Skerma_CRL_canola = lme(LogSkermanella~ LogCoarseRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Skerma_TRL_intc) # 691.7127
AIC(Skerma_CRL_min) #  688.2221
AIC(Skerma_CRL_full)  # 689.6539 
AIC(Skerma_CRL_canola)# 690.2221

summary(Skerma_TRL_intc)
summary(Skerma_CRL_min)
summary(Skerma_CRL_full)
summary(Skerma_CRL_canola)

# LogStenotrophomonas

Stenot_TRL_intc = lm(LogStenotrophomonas ~ 1, data=vmydata)
Stenot_TRL_min =  lm(LogStenotrophomonas~ LogTotal_RL, data=vmydata)
Stenot_TRL_full = lme(LogStenotrophomonas~ LogTotal_RL , random=~1|week, method ="ML", data=vmydata)
Stenot_TRL_canola = lme(LogStenotrophomonas~ LogTotal_RL , random=~1|cline, method ="ML", data=vmydata)

AIC(Stenot_TRL_intc) # 864.2343
AIC(Stenot_TRL_min) # 841.5513
AIC(Stenot_TRL_full) # 807.3586
AIC(Stenot_TRL_canola)# 843.5425

summary(Stenot_TRL_intc)
summary(Stenot_TRL_min)
summary(Stenot_TRL_full)
summary(Stenot_TRL_canola)

#
Stenot_TRL_intc = lm(LogStenotrophomonas ~ 1, data=vmydata)
Stenot_TFRL_min =  lm(LogStenotrophomonas~ LogTFRL2mm, data=vmydata)
Stenot_TFRL_full = lme(LogStenotrophomonas~ LogTFRL2mm , random=~1|week, method ="ML", data=vmydata)
Stenot_TFRL_canola = lme(LogStenotrophomonas~ LogTFRL2mm , random=~1|cline, method ="ML", data=vmydata)

AIC(Stenot_TRL_intc) # 864.2343
AIC(Stenot_TFRL_min) # 843.5201
AIC(Stenot_TFRL_full) # 807.3094   
AIC(Stenot_TFRL_canola)# 845.5175

summary(Stenot_TRL_intc)
summary(Stenot_TFRL_min)
summary(Stenot_TFRL_full)
summary(Stenot_TFRL_canola)

#

Stenot_TRL_intc = lm(LogStenotrophomonas ~ 1, data=vmydata)
Stenot_FRL_min =  lm(LogStenotrophomonas~ LogFRL, data=vmydata)
Stenot_FRL_full = lme(LogStenotrophomonas~ LogFRL , random=~1|week, method ="ML", data=vmydata)
Stenot_FRL_canola = lme(LogStenotrophomonas~ LogFRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Stenot_TRL_intc) # 864.2343
AIC(Stenot_FRL_min) # 838.2514
AIC(Stenot_FRL_full) # 807.5404    
AIC(Stenot_FRL_canola)# 840.1386

summary(Stenot_TRL_intc)
summary(Stenot_FRL_min)
summary(Stenot_FRL_full)
summary(Stenot_FRL_canola)

#
Stenot_TRL_intc = lm(LogStenotrophomonas ~ 1, data=vmydata)
Stenot_EFRL_min =  lm(LogStenotrophomonas~ LogExtraFRL, data=vmydata)
Stenot_EFRL_full = lme(LogStenotrophomonas~ LogExtraFRL , random=~1|week, method ="ML", data=vmydata)
Stenot_EFRL_canola = lme(LogStenotrophomonas~ LogExtraFRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Stenot_TRL_intc) # 864.2343
AIC(Stenot_EFRL_min) # 847.3551
AIC(Stenot_EFRL_full) # 807.192   
AIC(Stenot_EFRL_canola)# 849.3551

summary(Stenot_TRL_intc)
summary(Stenot_EFRL_min)
summary(Stenot_EFRL_full)
summary(Stenot_EFRL_canola)

#

Stenot_TRL_intc = lm(LogStenotrophomonas ~ 1, data=vmydata)
Stenot_CRL_min =  lm(LogStenotrophomonas~ LogCoarseRL, data=vmydata)
Stenot_CRL_full = lme(LogStenotrophomonas~ LogCoarseRL , random=~1|week, method ="ML", data=vmydata)
Stenot_CRL_canola = lme(LogStenotrophomonas~ LogCoarseRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Stenot_TRL_intc) # 864.2343
AIC(Stenot_CRL_min) #  830.7127
AIC(Stenot_CRL_full)  # 689.6539 
AIC(Stenot_CRL_canola) # 832.7127

summary(Stenot_TRL_intc)
summary(Stenot_CRL_min)
summary(Stenot_CRL_full)
summary(Stenot_CRL_canola)


# LogVariovorax

Variov_TRL_intc = lm(LogVariovorax ~ 1, data=vmydata)
Variov_TRL_min =  lm(LogVariovorax~ LogTotal_RL, data=vmydata)
Variov_TRL_full = lme(LogVariovorax~ LogTotal_RL , random=~1|week, method ="ML", data=vmydata)
Variov_TRL_canola = lme(LogVariovorax~ LogTotal_RL , random=~1|cline, method ="ML", data=vmydata)

AIC(Variov_TRL_intc) # 831.2716
AIC(Variov_TRL_min) # 797.3162
AIC(Variov_TRL_full) # 776.3505
AIC(Variov_TRL_canola)# 799.13

summary(Variov_TRL_intc)
summary(Variov_TRL_min)
summary(Variov_TRL_full)
summary(Variov_TRL_canola)

#
Variov_TRL_intc = lm(LogVariovorax ~ 1, data=vmydata)
Variov_TFRL_min =  lm(LogVariovorax~ LogTFRL2mm, data=vmydata)
Variov_TFRL_full = lme(LogVariovorax~ LogTFRL2mm , random=~1|week, method ="ML", data=vmydata)
Variov_TFRL_canola = lme(LogVariovorax~ LogTFRL2mm , random=~1|cline, method ="ML", data=vmydata)

AIC(Variov_TRL_intc) # 831.2716
AIC(Variov_TFRL_min) # 798.7099
AIC(Variov_TFRL_full) # 776.4711   
AIC(Variov_TFRL_canola)# 800.555

summary(Variov_TRL_intc)
summary(Variov_TFRL_min)
summary(Variov_TFRL_full)
summary(Variov_TFRL_canola)

#

Variov_TRL_intc = lm(LogVariovorax ~ 1, data=vmydata)
Variov_FRL_min =  lm(LogVariovorax~ LogFRL, data=vmydata)
Variov_FRL_full = lme(LogVariovorax~ LogFRL , random=~1|week, method ="ML", data=vmydata)
Variov_FRL_canola = lme(LogVariovorax~ LogFRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Variov_TRL_intc) # 831.2716
AIC(Variov_FRL_min) # 794.0276
AIC(Variov_FRL_full) # 775.7645    
AIC(Variov_FRL_canola)# 795.0428

summary(Variov_TRL_intc)
summary(Variov_FRL_min)
summary(Variov_FRL_full)
summary(Variov_FRL_canola)

#
Variov_TRL_intc = lm(LogVariovorax ~ 1, data=vmydata)
Variov_EFRL_min =  lm(LogVariovorax~ LogExtraFRL, data=vmydata)
Variov_EFRL_full = lme(LogVariovorax~ LogExtraFRL , random=~1|week, method ="ML", data=vmydata)
Variov_EFRL_canola = lme(LogVariovorax~ LogExtraFRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Variov_TRL_intc) # 831.2716
AIC(Variov_EFRL_min) # 803.8701
AIC(Variov_EFRL_full) # 777.0659   
AIC(Variov_EFRL_canola)# 805.8608

summary(Variov_TRL_intc)
summary(Variov_EFRL_min)
summary(Variov_EFRL_full)
summary(Variov_EFRL_canola)

#

Variov_TRL_intc = lm(LogVariovorax ~ 1, data=vmydata)
Variov_CRL_min =  lm(LogVariovorax~ LogCoarseRL, data=vmydata)
Variov_CRL_full = lme(LogVariovorax~ LogCoarseRL , random=~1|week, method ="ML", data=vmydata)
Variov_CRL_canola = lme(LogVariovorax~ LogCoarseRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Variov_TRL_intc) # 831.2716
AIC(Variov_CRL_min) #  797.9092
AIC(Variov_CRL_full)  # 776.808 
AIC(Variov_CRL_canola)# 799.6202

summary(Variov_TRL_intc)
summary(Variov_CRL_min)
summary(Variov_CRL_full)
summary(Variov_CRL_canola)

# LogXanthomonas

Xantho_TRL_intc = lm(LogXanthomonas ~ 1, data=vmydata)
Xantho_TRL_min =  lm(LogXanthomonas~ LogTotal_RL, data=vmydata)
Xantho_TRL_full = lme(LogXanthomonas~ LogTotal_RL , random=~1|week, method ="ML", data=vmydata)
Xantho_TRL_canola = lme(LogXanthomonas~ LogTotal_RL , random=~1|cline, method ="ML", data=vmydata)

AIC(Xantho_TRL_intc) # 293.6256
AIC(Xantho_TRL_min) # 295.6206
AIC(Xantho_TRL_full) # 297.6206
AIC(Xantho_TRL_canola)# 297.6206

summary(Xantho_TRL_intc)
summary(Xantho_TRL_min)
summary(Xantho_TRL_full)
summary(Xantho_TRL_canola)

#
Xantho_TRL_intc = lm(LogXanthomonas ~ 1, data=vmydata)
Xantho_TFRL_min =  lm(LogXanthomonas~ LogTFRL2mm, data=vmydata)
Xantho_TFRL_full = lme(LogXanthomonas~ LogTFRL2mm , random=~1|week, method ="ML", data=vmydata)
Xantho_TFRL_canola = lme(LogXanthomonas~ LogTFRL2mm , random=~1|cline, method ="ML", data=vmydata)

AIC(Xantho_TRL_intc) # 293.6256
AIC(Xantho_TFRL_min) # 295.6091
AIC(Xantho_TFRL_full) # 297.6091   
AIC(Xantho_TFRL_canola)# 297.6091

summary(Xantho_TRL_intc)
summary(Xantho_TFRL_min)
summary(Xantho_TFRL_full)
summary(Xantho_TFRL_canola)

#

Xantho_TRL_intc = lm(LogXanthomonas ~ 1, data=vmydata)
Xantho_FRL_min =  lm(LogXanthomonas~ LogFRL, data=vmydata)
Xantho_FRL_full = lme(LogXanthomonas~ LogFRL , random=~1|week, method ="ML", data=vmydata)
Xantho_FRL_canola = lme(LogXanthomonas~ LogFRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Xantho_TRL_intc) # 293.6256
AIC(Xantho_FRL_min) # 295.6091
AIC(Xantho_FRL_full) # 297.6091    
AIC(Xantho_FRL_canola)# 297.6091

summary(Xantho_TRL_intc)
summary(Xantho_FRL_min)
summary(Xantho_FRL_full)
summary(Xantho_FRL_canola)

#
Xantho_TRL_intc = lm(LogXanthomonas ~ 1, data=vmydata)
Xantho_EFRL_min =  lm(LogXanthomonas~ LogExtraFRL, data=vmydata)
Xantho_EFRL_full = lme(LogXanthomonas~ LogExtraFRL , random=~1|week, method ="ML", data=vmydata)
Xantho_EFRL_canola = lme(LogXanthomonas~ LogExtraFRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Xantho_TRL_intc) # 293.6256
AIC(Xantho_EFRL_min) # 295.5849
AIC(Xantho_EFRL_full) # 297.5849   
AIC(Xantho_EFRL_canola)# 297.5849

summary(Xantho_TRL_intc)
summary(Xantho_EFRL_min)
summary(Xantho_EFRL_full)
summary(Xantho_EFRL_canola)

#

Xantho_TRL_intc = lm(LogXanthomonas ~ 1, data=vmydata)
Xantho_CRL_min =  lm(LogXanthomonas~ LogCoarseRL, data=vmydata)
Xantho_CRL_full = lme(LogXanthomonas~ LogCoarseRL , random=~1|week, method ="ML", data=vmydata)
Xantho_CRL_canola = lme(LogXanthomonas~ LogCoarseRL , random=~1|cline, method ="ML", data=vmydata)

AIC(Xantho_TRL_intc) # 293.6256
AIC(Xantho_CRL_min) #  295.4126
AIC(Xantho_CRL_full)  # 297.4126 
AIC(Xantho_CRL_canola)# 297.4126

summary(Xantho_TRL_intc)
summary(Xantho_CRL_min)
summary(Xantho_CRL_full)
summary(Xantho_CRL_canola)





