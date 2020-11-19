#BENEFICIAL PLANT??MICROBE INTERACTION IN AGROECOSYSTEMS: DECIPHERING THE RHIZOSPHERE MICROBIAL COMMUNITY 
#IN FIELD GROWN Brassica napus L.

# 1.	Sequence data processing

## FASTQ files were processed using DADA2 implemented in QIIME2. Data processing was done on a server.
# Below you will find script to install and get started with your sequence data processiong.

# a.	Set up for processing sequences reads using QIIME2 on a server  

## First obtain access to server
##ssh to the server and provide your password as needed

ssh xxx@oats.usask.ca # ssh to your account name xxx, provide password when prompted.

## create working director for your analysis. Remember where you created your directory.

mkdir canola_2016

## Install Miniconda and qiime2: using the guide for Linux system. Please refer to QIIME2 website 
## (https://qiime2.org/) for more details.

## install miniconda first

wget https://repo.continuum.io/miniconda/Miniconda3-3.7.0-Linux-x86_64.sh -O ~/miniconda.sh

bash ~/miniconda.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH"
conda update conda
conda install wget

## install qiime2: through the Linux approach from qiime2 webpage. See QIIME2 webpage 
## (qiime2.https://org/) for more details if needed.

wget https://data.qiime2.org/distro/core/qiime2-2019.1-py36-linux-conda.yml
conda env create -n qiime2-2019.1 --file qiime2-2019.1-py36-linux-conda.yml

## Useful scripts for working between your local machine and server
## Transferring files to and from server. First exist from the server and cd to your local directory 
## that contains files you want to transfer.

## to transfer individual files

scp paired-end-trimmed_MF_17.qza xxx@oats:/home/xxx/canola_2016
scp manifest.csv xxx@oats:/home/xxx/canola_2016
scp MF_metadata.txt xxx@oats:/home/xxx/canola_2016

## to transfer the whole folder

scp -r Soil xxx@oats:/home/xxx/canola_2016/soilraw
scp -r primer_trimmed_fastqs xxx@oats:/home/xxx/canola_2016

## To transfer from server to local folder

scp xxx@oats:/home/xxx/canola_2016/paired-end-trimmed.qza /mnt/d
scp xxx@oats:/home/xxx/canola_2016/paired-end-trimmed.qzv /mnt/d

## to copy a file from one directory into another in linux

cp my_SILVA_V3_V5_342F_806R_qiime2_2019_1_classifier.qza /home/xxx/canola_2016/
  my_SILVA_V3_V5_342F_806R_qiime2_2019_1_classifier.qza

## to extract tar.gz files

tar xvzf gg_13_8_otus.tar.gz

#### ................................#### ........................................... ####

## b.	Process the sequence data using QIIME2

## import your data (fastq) into qiime. To import an already demultiplexed sequence files use the manifest 
## approach in qiime tutorial (https://qiime2.org/).

## use qiime tools import command to import. Copy and paste the following command. J
## ust change the input-path and output- path for your specific case

qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path manifest.csv \
--output-path paired-end-trimmed.qza \
--input-format PairedEndFastqManifestPhred33

## visualized the imported paired end reads
qiime demux summarize \
--i-data paired-end-trimmed.qza \
--o-visualization paired-end-trimmed.qzv    ## qzv files can be visualized in qiime2's webpage (qiime2.https://org/)

## denoising DAAD2: trim and trunc specification will be decided based on quality score. 
## Examine the qza file and decide on these parameters. You will also be able to 
## check how much of your sequence were retained after applying your parameters (denoising stats).

qiime dada2 denoise-paired \
--i-demultiplexed-seqs paired-end-trimmed.qza \
--p-trim-left-f 13 \
--p-trim-left-r 0 \
--p-trunc-len-f 260 \
--p-trunc-len-r 208 \
--o-table table.qza \
--o-representative-sequences rep-seqs.qza \
--o-denoising-stats denoising-stats.qza

## visualize the denoising stats

qiime metadata tabulate \
--m-input-file denoising-stats.qza \
--o-visualization denoising-stats.qzv

##FeatureTable and FeatureData summaries: generate summaries of the artifacts generate by dada2 denoising step

qiime feature-table summarize \
--i-table table.qza \
--o-visualization table.qzv \
--m-sample-metadata-file metadata.txt

qiime tools export \
--input-path table.qza \
--output-path feature-table
qiime feature-table tabulate-seqs \
--i-data rep-seqs.qza \
--o-visualization rep-seqs.qzv

## Now train taxonomy classifier and assign taxonomy
## train classifying sequence.
## Download latest greengenes database 
## use the 99 seq and taxonomy for training
## change sequence in to qiime artifact

qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path 99_otus.fasta \
--output-path 99_seqs.qza 

qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path 99_otu_taxonomy.txt \
--output-path ref-taxonomy.qza

qiime feature-classifier extract-reads \
--i-sequences 99_seqs.qza \
--p-f-primer CTACGGGGGGCAGCAG  \
--p-r-primer GGACTACCGGGGTATCT \
--p-min-length 100 \
--p-max-length 440 \
--o-reads ref-seqsg.qza

qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads ref-seqsg.qza \
--i-reference-taxonomy ref-taxonomy.qza \
--o-classifier my_ggenes_V3_V5_342F_806R_qiime2_2019_1_classifier.qza

## Now you have trained your classifier proceed to testing it on your denoised rep.sequence.

qiime feature-classifier classify-sklearn \
--i-classifier my_ggenes_V3_V5_342F_806R_qiime2_2019_1_classifier.qza \
--i-reads rep-seqs.qza \
--o-classification taxonomy2g.qza

qiime tools export \
--input-path taxonomy2g.qza \
--output-path taxonomy-with-spaces2g
qiime metadata tabulate \
--m-input-file taxonomy-with-spaces2g/taxonomy.tsv  \
--o-visualization taxonomy-as-metadata2.qzv

qiime tools export \
--input-path taxonomy-as-metadata2.qzv \
--output-path taxonomy-as-metadata2g

qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-path taxonomy-as-metadata2g/metadata.tsv \
--output-path taxonomy-without-spaces2g.qza

qiime metadata tabulate \
--m-input-file taxonomy-without-spaces2g.qza \
--o-visualization taxonomy-without-spaces2g.qzv

qiime tools export \
--input-path taxonomy-without-spaces2g.qza \
--output-path taxonomy-without-spaces2g

## generate interactive bar plots with the taxonomy

qiime taxa barplot \
--i-table table.qza \
--i-taxonomy taxonomy-without-spaces2g.qza \
--m-metadata-file metadata.txt \
--o-visualization taxa-bar-plots2g.qzv

## Filtering mitochondria, chloroplast from: table.qza, rep.seq
#filter table-qza

qiime taxa filter-table \
--i-table table.qza \
--i-taxonomy taxonomy-without-spaces2g.qza \
--p-exclude mitochondria,chloroplast \
--o-filtered-table table-no-mitochondria-no-chloroplast2g.qza

qiime feature-table summarize \
--i-table table-no-mitochondria-no-chloroplast2g.qza \
--o-visualization table-no-mitochondria-no-chloroplast2g.qzv \
--m-sample-metadata-file metadata.txt

## export filtered feature table as biome
qiime tools export \
--input-path table-no-mitochondria-no-chloroplast2g.qza \
--output-path feature-table-m-c-filtered2g

##Generate a tree 

qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences rep-seqs.qza \
--o-alignment aligned-rep-seqs.qza \
--o-masked-alignment masked-aligned-rep-seqs.qza \
--o-tree unrooted-tree.qza \
--o-rooted-tree rooted-tree.qza











