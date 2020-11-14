#BENEFICIAL PLANT??MICROBE INTERACTION IN AGROECOSYSTEMS: DECIPHERING THE RHIZOSPHERE MICROBIAL COMMUNITY 
#IN FIELD GROWN Brassica napus L.

# c.	QIIME2 to phyloseq object: example script to move from Qiime2 output to working with phyloseq object

setwd("D:/.. ") # set working directory

## read the otu table. This is qiime2 feature tabel (asv) converted to txt format. Conver the biome to txt file.

##Mitochondria and chloroplast are removed 

otu<-read.table(file = "otu_table-no-mc_SC.txt",header = TRUE)
head(otu)

##read taxonomy ## from qiime2

tax<-read.table(file = "taxonomy_SC.txt",sep = '\t', header = TRUE)
head(tax)

## merge the above two files: the merging is necessary as the two files have different row length. 
## We have filtered out chloroplast and mitochondria from the asv but not from taxonomy. By merging we will match them.

merged_file<-merge(otu,tax,by.x=c("OTUID"),by.y=c("OTUID"))

## output merged file

write.table(merged_file, file = "combined_otu_tax_SC",sep = '\t',col.names = TRUE,row.names = FALSE)

## You need to open the merged files and create two files. 
## one for taxonomy taking the seq and taxonomy columns (also change the ';' separated taxon to column).
## This is done manually

## asv table with sample names and abundance of each sequence. 
## Then these independent files will be imported to R again for further processing.

### upload the final tax and seq(otu) files and convert into phyloseq object

## load the following libraries or install if you don't have

library(ggplot2)
library(phyloseq)
library(ape)
library(biomformat)

## read seq table

otu_table=read.csv("seq_final_SC.csv", sep = ",",row.names = 1)
otu_table=as.matrix(otu_table)

## read in taxonomy: have been separated into different taxon level above- manually in excel

taxonomy=read.csv("tax_final_SC.csv",sep="," ,row.names=1)
taxonomy=as.matrix(taxonomy)

## metadata=read.table("metadata.txt",row.names = 1)# Reading it as txt file created a problem downstream 
## by not showing sample variables in phyloseq object

metadata=read.csv("SC_metadata.csv",sep="," ,row.names=1)
head(metadata)

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

physeq=phyloseq(OTU,TAX,META)

sample_variables(physeq)## check if sample variables are correctly presented in the phyloseq object

## now collect the Aliivibrio and adjust for sequencing depth

can_rhi_2016.fischeri <- subset_taxa(physeq, Genus == "Aliivibrio" )

## remove internal standard from dataset

can_rhi_2016.no_fischeri <- subset_taxa(physeq, Genus != "Aliivibrio")

## now using the internal standard to correct the abundance

fischeri <- as.vector(sample_sums(can_rhi_2016.fischeri))   ## vector created for the internal standard
wi <- 1e-10  ## weight of internal standard gDNA added to the samples
gi <- 4.49e-15  ### Weight of genome of the internal standard
ci <- 8   ##16S copy number of the internal standard
adj <- wi/gi*ci
fischeri.1 <- adj/fischeri
fischeri.1[is.infinite(fischeri.1)] <- 0
snv.table <- as.matrix(can_rhi_2016.no_fischeri@otu_table)  ### export asv table as matrix
snv.table.norm <- as.data.frame(sweep(snv.table, 2, fischeri.1, `*`))  # This divides the abundance matrix 
                                                                      # by the vector to normalize the data
tax.table.norm <- as.matrix(can_rhi_2016.no_fischeri@tax_table)  # This will extract taxonomic information 
                                                                # excluding the internal standard 

## recreate the phyloseq object using the internal standard abundance adjusted otutable
## import them as phyloseq objects

OTU=otu_table(snv.table.norm,taxa_are_rows=TRUE)
TAX=tax_table(tax.table.norm)
META=sample_data(metadata)

## check that your otu names are consistent across objects

taxa_names(TAX)
taxa_names(OTU)

## make sure files have the same sample names

sample_names(OTU)
sample_names(META)

## Everything looks good --- now proceed to merging them into one phyloseq object

physeq.norm=phyloseq(OTU,TAX,META)
sample_variables(physeq.norm)

## Further processing may be needed based on your data. For example to remove duplicates (PCR and extraction), checks, etc.
## You can subset by week separate fiels or you can during analysis within your phyloseq object 

rhiz_2016.norm.w1 <- subset_samples(can_2016.no.dup.2, Week == "1")   # example Subset to soil w1 samples

rhiz_2016.norm.2.w1<- prune_taxa(taxa_sums(rhiz_2016.norm.w1) > 0, rhiz_2016.norm.w1) # Remove zero-sum taxa before                                                                                             # outputting

## Remove zero-sum taxa before outputting (this also makes subsequent processing steps proceed much more quickly)

can_2016.norm.soil_zeor_fintered <- prune_taxa(taxa_sums(can_2016.no.dup.2) > 0, can_2016.no.dup.2)

## save as phyloseq object the zero sum filtered object
save(can_2016.norm.soil_zeor_fintered,file="D:/out_puts/..../myphylo_zero_sum_filtered_SC_17.RData")

## export taxa, asv and metadata files of the filtered phyloseq object

can.snv.table.2 <- as.matrix(can_2016.norm.soil_zeor_fintered@otu_table)
can.tax.table.2 <- as.matrix(can_2016.norm.soil_zeor_fintered@tax_table)
can.sam.table.2 <- as.data.frame(can_2016.norm.soil_zeor_fintered@sam_data)

# write the above files
write.csv(can.snv.table.2, "D:/out_puts/...../can.asv.rhizo_SC.csv")
write.csv(can.tax.table.2, "D:/out_puts/...../can.tax.rhizo_SC.csv")
write.csv(can.sam.table.2, "D:/out_puts/....../can.sam.rhizo_SC.csv")

##creat biomfile: for future use when needed

can.2016.biom.rhizo <- make_biom(can.snv.table.2, can.sam.table.2, can.tax.table.2)
write_biom(can.2016.biom.rhizo, "D:/out_puts//canola_rhizo_SC.biom")
















