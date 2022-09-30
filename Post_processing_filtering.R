# TITLE: FILTERING FOR SNPs ANALYSES
# Authors: Bernat Burriel-Carranza & Loukia Spilani
# Date: 20/08/2020

# Set the working directory to where you have the vcf ipyrad files.
setwd("PATH_TO_FILE")

# Load packages
library("radiator");library("dplyr");library("gdsfmt");library("SNPRelate");library("vcfR");library("tidyverse");library("SeqArray")


# 1. Prepare a strata file
# It’s a tab separated file, e.g. example.strata.tsv. This file has to be created outside R
# It has to be prepared with the species codes of each individual present in the ipyrad dataset and it a number assigning each individual to a strata (population or group) can be given
# A minimum of 2 columns: INDIVIDUALS and STRATA is required.
# The STRATA column identifies the individuals stratification, the hierarchical groupings: populations, sampling sites or any grouping you want.
# It’s like stacks population map file with header…
## Load your .vcf file and blacklist the individuals that are not desired
  
## Read and inspect the STRATA-FILE
strata <- read.table('strata_file.tsv',sep = "\t", header = TRUE)
head(strata)
tail(strata)
# If you only want specific populations subset the populations object
populations <- unique(strata[,2])
populations

## Read the vcf file and select the desired strata. Here we do not further filter, we will do it in the following steps
data <- read_vcf("<.path/to/vcf>",
                 strata = "<path/to/strata_file>",
                 interactive.filter = FALSE,
                 path.folder = 'my_folder', 
                 filter.common.markers = FALSE,
                 filter.monomorphic = FALSE)


#pop.select = populations["desired_populaitons"], # add this chunk above to select populations

## Filter for those markers that are present in at least one individual of each strata_group
data <- filter_common_markers(data = data,
                             path.folder = 'my_folder')


## Decide the missigness that you want to allow
missingness <- c(0.98, 0.96, 0.94, 0.92,0.90)#,0.88,0.86,0.84,0.82,0.80,0.78,0.76,0.74,0.72,0.70, 0.68,0.66,0.64,0.62,0.60)
#                 0.58,0.56,0.54,0.52,0.50, 0.48,0.46,0.44,0.42,0.40)#0.38,0.36,0.34,0.32, 0.30)#, 0.28,0.26, 0.24,0.22 ,0.20,0.18, 0.16,0.14,0.12,0.1)

## Create   a loop where iteratively, individual and genotiping filterings are aplied up to the desired percentage of missigness.
data1 <- data
for (i in 1:length(missingness)) {
  data1 <- filter_individuals(data = data1,
                              interactive.filter = FALSE,
                              filter.individuals.missing = missingness[i],
                              filter.individuals.heterozygosity = NULL,
                              filter.individuals.coverage.total = NULL,
                              path.folder = 'my_folder')
  data1 <- filter_genotyping(data = data1,
                             interactive.filter = FALSE,
                             filter.genotyping = missingness[i],
                             path.folder = 'my_folder')
}

# Apply a hard genotyping filter if you want less missing data in your final dataset while retaining more individuals
data1 <- filter_genotyping(data = data1,
                           interactive.filter = FALSE,
                           filter.genotyping = 0.2,
                           path.folder = 'my_folder')


## Apply a minor allele count filtering. The Minor Allele Count function remove/blacklist markers based on Minor/Alternate Allele Count (MAC).
#   Use it to remove noise, sequencing errors or low polymorphism markers. Some analysis performs better with the full spectrum of allele frequency,
#   so careful with high threshold that inevitably results in biaises. Leave it at 1
data1 <- filter_mac(data = data1,
                    interactive.filter = FALSE,
                    path.folder = 'my_folder',
                    filter.mac = 1)
# Apply a minot allel frequency filter: Minor Allele Filter. Remove markers based on Minor/Alternate Allele Frequency (MAF) or Count (MAC). VALUE = 0.05
data1 <- filter_maf(data = data1,
                    interactive.filter = FALSE,
                    path.folder = 'my_folder',
                    maf.thresholds = c('SNP', NULL, 'OR', 0.05, 1))

## Remove monomorphic, non-informative-sites
data1 <- filter_monomorphic(data = data1,
                            path.folder = 'my_folder',
                            verbose = TRUE)

# tidy your data to export it as a VCF
tidy_data <- tidy_genomic_data(data = data1)
position <- data.frame(MARKERS = tidy_data$MARKERS) %>%
  separate(col = MARKERS, into = c('Loci', 'numloc', 'space', 'pos', 'space2', 'numpos'), sep = '_')
position$Loci <- 'loc'
position$pos <- 'pos'
ids <- paste(position$Loci, position$numloc, '_', position$pos, position$numpos, sep = '')
tidy_data$LOCUS <- ids

# remove biallellic sites
tidy_data_biallellic <- tidy_data[grep(',', tidy_data$ALT, invert = TRUE),] 
tidy_data_biallellic <- tidy_data_biallellic[grep(',', tidy_data_biallellic$REF, invert = TRUE),] 
write_vcf(data = tidy_data_biallellic,
          filename = 'my_folder/ind_gen_filetering_dataset')

system("plink2 --vcf my_folder/ind_gen_filetering_dataset.vcf --const-fid --allow-extra-chr --max-alleles 2 -min-alleles 2 --make-bed --out my_folder/concatenated_snps_biallellic_only") #bi=bi-allelic only
system("plink2 --allow-extra-chr --bfile my_folder/concatenated_snps_biallellic_only --export vcf-4.2 --out my_folder/FILE_NAME") # Export into vcf to obtain the Concatenated SNPs

#########################################################################################################################################################
#### 2. Create a SNP list to filter out those SNPs from a locus that are present in less samples (NS) based on the remaining snps after the
#### filtering-------------------------------------------------------------------------

## Create a gds file from the Plink output to work with SNPRelate
##load the plink files to R
bed1.fn <- "my_folder/concatenated_snps_biallellic_only.bed"
fam1.fn <- "my_folder/concatenated_snps_biallellic_only.fam"
bim1.fn <- "my_folder/concatenated_snps_biallellic_only.bim"

##cvt.chr: "int" - chromosome code in the GDS file is integer; "char" - chromosome code in the GDS file is character
SNPRelate::snpgdsBED2GDS(bed1.fn, fam1.fn, bim1.fn, cvt.chr="char", "my_folder/concatenated_snps_biallellic_only.gds")

## Open the GDS file
genofile1 <- SNPRelate::snpgdsOpen("my_folder/concatenated_snps_biallellic_only.gds")
#snpgdsSummary(genofile1)
chr1 <- read.gdsn(index.gdsn(genofile1, "snp.chromosome"))
pos1 <- read.gdsn(index.gdsn(genofile1, "snp.position"))
id1 <- read.gdsn(index.gdsn(genofile1, "snp.id"))
id1.v <- as.vector(unlist(unname(id1)))

## Import VCF
my.vcf <- vcfR::read.vcfR('my_folder/ind_gen_filetering_dataset.vcf')

##Read all the fixed columns of a vcf (CHROM, POS etc. as a data frame)
##Both chromR objects and vcfR objects contain a region with fixed variables. getFIX allows you to isolate these variables
## from these objects while virtually ignoring the actual genetic data. INFO2df reformats INFO data as a data.frame and
##  handles class when possible.
my.vcf.df <- cbind(as.data.frame(vcfR::getFIX(my.vcf)), vcfR::INFO2df(my.vcf))
my.vcf.df <- my.vcf.df[grep(',', my.vcf.df$ALT, invert = TRUE),] 
my.vcf.df <- my.vcf.df[grep(',', my.vcf.df$REF, invert = TRUE),] 

##Filter the SNP list based on the remaining SNP IDs after the missingnes and biallelic
filt.vcf.df <- filter(my.vcf.df, ID %in% id1.v)

##Group by locus (CHROM). From each group (=locus) keep the rows (=SNPs) that have max number of samples (NS).
##You get all the SNPs per locus that have max NS. If all SNPs on a locus have equal number of NS it keeps all of them
filt.vcf.df <- dplyr::group_by(filt.vcf.df, CHROM, .add = FALSE) %>%
  filter(NS == max(NS)) %>%
  dplyr::arrange(CHROM, POS)

##Ungroup the data frame
filt.vcf.df <- dplyr::ungroup(filt.vcf.df, CHROM)

##Extract SNP IDs column to a txt
snpset <- filt.vcf.df %>% dplyr::select(ID)
snpset.id <- unlist(unname(snpset))
write.table(snpset.id, "my_folder/snpset.id.txt", row.names = F, col.names = F, quote = FALSE, sep="\t")

## We use 'plink' to extract the set of slected SNPs based on Number of Samples (see section above)
system("plink2 --allow-extra-chr --bfile my_folder/concatenated_snps_biallellic_only --extract my_folder/snpset.id.txt --max-alleles 2 --make-bed --out my_folder/test_bi_cr_maf_snpNS") #snpNS: snps with less number of samples per locus

####3. Select one snps per locus at random------------------------------------------------
## Create a gds file from the Plink output to work with SNPRelate
##load the plink files to R
bed2.fn <- "my_folder/test_bi_cr_maf_snpNS.bed"
fam2.fn <- "my_folder/test_bi_cr_maf_snpNS.fam"
bim2.fn <- "my_folder/test_bi_cr_maf_snpNS.bim"

##cvt.chr: "int" - chromosome code in the GDS file is integer; "char" - chromosome code in the GDS file is character
SNPRelate::snpgdsBED2GDS(bed2.fn, fam2.fn, bim2.fn, cvt.chr="char", "my_folder/test_bi_cr_maf_snpNS.gds")

## Open the GDS file
genofile2 <- SNPRelate::snpgdsOpen("my_folder/test_bi_cr_maf_snpNS.gds")
#snpgdsSummary(genofile2)

chr2 <- read.gdsn(index.gdsn(genofile2, "snp.chromosome"))
pos2 <- read.gdsn(index.gdsn(genofile2, "snp.position"))
id2 <- read.gdsn(index.gdsn(genofile2, "snp.id"))
id2.v <- as.vector(unlist(unname(id2)))

#set seed so that I can replicate the anaysis
#select one snp per chromosome randomly (gives back a list of logical values)
set.seed(447)
randsnpset <- snpgdsApartSelection(chr2, pos2, max.n.snp.perchr=1, verbose=FALSE)

## transform the list of logical values to a dtaframe and saving it for cross checking
randsnpset.logic <- as.data.frame(unlist(unname(randsnpset)))
write.table(randsnpset.logic, "my_folder/randsnpset.logic.txt", row.names = F, col.names = F, quote = FALSE, sep="\t")

##based on the list of logical values extract from the SNP id list those snps for which the logical value is TRUE
##write and save the table of selected snp id so that they can be extracted from the vcf file
selectedIDs <- id2.v[randsnpset]
write.table(selectedIDs, "my_folder/randSelectedIDs.txt", row.names = F, col.names = F, quote = FALSE, sep="\t")

## We use 'plink' to extract the set of randomly slected SNPs
system("plink2 --allow-extra-chr --bfile my_folder/test_bi_cr_maf_snpNS --max-alleles 2 --extract my_folder/randSelectedIDs.txt --make-bed --out my_folder/qc_FILE_NAME") #randsnp=random snp

## Export it in 012 format
system("plink --bfile my_folder/qc_FILE_NAME --recode12 --out my_folder/qc_FILE_NAME_12 --allow-extra-chr") #randsnp=random snp

## Close the genotype file
snpgdsClose(genofile1)
snpgdsClose(genofile2)

##Convert to .vcf
system("plink2 --allow-extra-chr --bfile my_folder/qc_FILE_NAME --export vcf-4.2 --out my_folder/qc_FILE_NAME_usnp")

rm(list = ls())
