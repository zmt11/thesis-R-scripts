### This scripts are to do several summary stats on the final phyloseq object
# load the phyloseq object saved as can.physq.data.2016
#### Summary statistics 2016 canola sequence data##
# load can.physq.data.2016 object

## Published article title: Core and Differentially Abundant Bacterial Taxa in the Rhizosphere of Field Grown Brassica #napus Genotypes: Implications for Canola Breeding

##2016 phyloseq object

can_2016.norm.soil_zeor_fintered

# There ate 12567 zero sum taxa in 477 samples
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 12567 taxa and 477 samples ]
#sample_data() Sample Data:       [ 477 samples by 22 sample variables ]
#tax_table()   Taxonomy Table:    [ 12567 taxa by 7 taxonomic ranks ]

#################### PART I #############################################

# number of taxa

ntaxa(can_2016.norm.soil_zeor_fintered) # there are 12567 taxa

#number of samples

nsamples(can_2016.norm.soil_zeor_fintered)## there are 477 samples 


### summery of the sample reads

sum(sample_sums(can_2016.norm.soil_zeor_fintered))## 1571433759 total number of reads in the dataset
min(sample_sums(can_2016.norm.soil_zeor_fintered))### 5925.293 minimum number of reads in the dataset
max(sample_sums(can_2016.norm.soil_zeor_fintered))### 354234808 maximum number of reads in the dataset
mean(sample_sums(can_2016.norm.soil_zeor_fintered))### 3294410 mean number of reads in the dataset

### creat a data table with total number of reads per sample (per canola line)

library(data.table)
sample.read.per.sample = data.table(as(sample_data(can_2016.norm.soil_zeor_fintered), "data.frame"),
                                    TotalReads = sample_sums(can_2016.norm.soil_zeor_fintered), keep.rownames = TRUE)

## write as csv fiel

write.csv(sample.read.per.sample, 'D:/out_puts/phyloseq/with_green_genes_taxa/summary_tables/totalread.per.sample_gg.csv')


## if you want to work on the rounded abundance valuece to nearest intiger run the following code first

canola2016_nearest_integer_gg<-ceiling(otu_table(can_2016.norm.soil_zeor_fintered))

### merge the otu table with integer values with taxa table and meta data

can.physq.data.2016.nearest.intiger_gg<-merge_phyloseq(canola2016_nearest_integer_gg,tax_table(can_2016.norm.soil_zeor_fintered),
                                                    sample_data(can_2016.norm.soil_zeor_fintered))

#save nearest intiger as RData for later easier loading

save(can.physq.data.2016.nearest.intiger_gg,file='D:/out_puts/phyloseq/with_green_genes_taxa/canola2016_nearest_integer_filtered_gg.RData')

## based on rounded to nearest intiger data

sum(sample_sums(can.physq.data.2016.nearest.intiger_gg))##1571475596 total number of reads in the dataset
min(sample_sums(can.physq.data.2016.nearest.intiger_gg))### 5936 minimum number of reads in the dataset
max(sample_sums(can.physq.data.2016.nearest.intiger_gg))### 354234976 maximum number of reads in the dataset
mean(sample_sums(can.physq.data.2016.nearest.intiger_gg))### 3294498 mean number of reads in the dataset

## the following is a function to get the most abundandant taxa for each sample and its corresponding abundance:find.top.taxa

find.top.taxa <- function(x,taxa){
  require(phyloseq)
  top.taxa <- tax_glom(x, taxa)
  otu <- otu_table(t(top.taxa))
  if (taxa_are_rows(otu)){
    otu <- t(otu)
  }
  tax <- tax_table(top.taxa)
  j<-apply(otu,1,which.max)
  k <- j[!duplicated(j)]
  l <- data.frame(tax[k,])
  m <- data.frame(otu[,k])
  colnames(m) = l[,taxa]
  n <- colnames(m)[apply(m,1,which.max)]
  m[,taxa] <- n
  return(m)
}
top.taxa <- find.top.taxa(can_2016.norm.soil_zeor_fintered,"Phylum")

#

### use the find.top.taxa function to get the most abundant phyla within each sample

top.taxa <- find.top.taxa(can_2016.norm.soil_zeor_fintered,"Phylum")

### save the top phyla per sample as csv file: summarize frequency of 

write.csv(top.taxa, 'D:/out_puts/phyloseq/with_green_genes_taxa/summary_tables/top.phylum.per.sample_gg.csv')


#### creat csv file with prevalence and total abunndance of each snv and add taxomomy of each snv

prevelancedf = apply(X = otu_table(can_2016.norm.soil_zeor_fintered),
                     MARGIN = 1,
                     FUN = function(x){sum(x > 0)})                              ### calculate prevalence of each snv


# Add taxonomy and total read counts to this data.frame

prevelancedf = data.frame(Prevalence = prevelancedf,
                          TotalAbundance = taxa_sums(can_2016.norm.soil_zeor_fintered),
                          tax_table(can_2016.norm.soil_zeor_fintered))
prevelancedf[1:10,]

# write it as a csv file

write.csv(prevelancedf, 'D:/out_puts/phyloseq/with_green_genes_taxa/summary_tables/Prevalenceandtotalabundace_each_snv_gg.csv')

## plot ##Plot of all phylum showing total abunandce and relative prevalence in each sample

library(ggplot2)

prevelancedf1 = subset(prevelancedf, Phylum %in% get_taxa_unique(can_2016.norm.soil_zeor_fintered, taxonomic.rank = "Phylum"))

ggplot(prevelancedf1, aes(TotalAbundance, Prevalence / nsamples(can_2016.norm.soil_zeor_fintered),color=Phylum)) +theme_bw()+
  geom_point(size = 2) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") + theme(axis.text.x=element_text(angle=90)) +
  facet_wrap(~Phylum) + theme(legend.position="none")

#prevelance/abundance phylum 

phylum_prevalece_totaladundance=plyr::ddply(prevelancedf, "Phylum", function(df1){
  data.frame(frequence=length(df1$Prevalence), mean_prevalence=mean(df1$Prevalence),total_abundance=sum(df1$TotalAbundance,na.rm = T),stringsAsFactors = F)
})

# write it as a csv file
write.csv(phylum_prevalece_totaladundance, 'D:/out_puts/phyloseq/with_green_genes_taxa/summary_tables/PHYLUM_Prevalence_totalabundace_gg.csv')


#### Summary of the three major phyla : total abundance, mean, Acidobacteria, Actinobacteria, Proteobacteria 

threemajorphyla = subset_taxa(can_2016.norm.soil_zeor_fintered, Phylum=="Acidobacteria" | Phylum=="Actinobacteria" 
                              |Phylum=="Proteobacteria")  ## subset the three phyla that you want to have summary for


threemajorphyla_glom=tax_glom(threemajorphyla, taxrank="Phylum")

threemajorphyla_glom_df = psmelt(threemajorphyla_glom)# change to data frame

summary_threemajorphyla=ddply(threemajorphyla_glom_df,c("CanolaLine","Phylum"), summarise, N=length(Abundance), 
                   Average_Abundace=mean(Abundance),Minimum=min(Abundance),Maximum=max(Abundance), Total=sum(Abundance))

### a data frame with the summary you want

write.csv(summary_threemajorphyla,"D:/out_puts/phyloseq/with_green_genes_taxa/summary_tables/summary_threemajorphyla_perline_gg.csv")# save this summary for reporting





