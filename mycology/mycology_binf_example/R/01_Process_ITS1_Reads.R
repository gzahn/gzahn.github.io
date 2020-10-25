# -----------------------------------------------------------------------------#
# Mediterranean sponge fungi - ocean acidification project 
# Processing cutadapt-trimmed Reads
# Author: Geoffrey Zahn
# Software versions:  ITSxpress v 1.8.0
#                     R v 3.6.3
#                     tidyverse v 1.3.0
#                     dada2 v 1.14.1
#                     decontam v 1.6.0
#                     phyloseq v 1.30.0
#                     purrr v 0.3.4
#                     Biostrings 2.54.0
# -----------------------------------------------------------------------------#

# EXTRACT ITS1 SEQUENCES ####

# Run ITSxpress on all forward fasta files (reverse reads aren't working well for some reason)
          
                                            # un-comment the following 2 lines if running this analysis from scratch!
# itsxpress <- "for i in ./seqs/fastqs/*R1_001.fastq.gz; do /home/gzahn/.local/bin/itsxpress --fastq $i  --outfile $i.FungalITS1.fastq.gz --region ITS1 --taxa Fungi -s --threads 4 --log $i.ITSx.log; done"
# system(command = itsxpress)
getwd()

# BiocManager::install("phyloseq",update = FALSE) # use this code to install phyloseq if problems arise

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(dada2); packageVersion("dada2")
library(decontam); packageVersion("decontam")
library(phyloseq); packageVersion("phyloseq")
library(purrr); packageVersion("purrr")
library(Biostrings); packageVersion("Biostrings")

source("./R/palettes.R")
source("./R/plot_bar2.R")


#################################################################################
#                               Main workflow                                   #
# Filter and trim, denoise, sample inferrence, chimera and contaminant removal, # 
# taxonomic assignment, combine sequence table and metadata                     #
#################################################################################

# PARSE FILE PATHS ####

# File parsing - For this, we will use only the forward illumina reads - make sure to move fwd reads into their own directory for simplest processing
path <- "./seqs/fastqs" # CHANGE to the directory containing your demultiplexed fastq files
filtpath <- file.path(path, "filtered") # Filtered files go into the filtered/ subdirectory
if(!file_test("-d", filtpath)) dir.create(filtpath) # make directory for filtered fqs if not already present
fns <- sort(list.files(file.path(path), full.names = TRUE, pattern = "_R1_001.fastq.gz.FungalITS1.fastq.gz"))


A <- unlist(map(strsplit(basename(fns), "\\."), 1))
B <- unlist(map(strsplit(basename(fns), "\\."), 2))
C <- unlist(map(strsplit(basename(fns), "\\."), 3))
sample.names <- paste(A,B,C,sep = "_") 
rm(list = c("A","B","C"))

# visualize a couple of fwd read quality profiles to help select reasonable filtration parameters
plotQualityProfile(fns[1:4])

# FILTER AND TRIM ####
filts <- file.path(path, "filtered", paste0(sample.names, "_filt.fastq.gz"))

out <- filterAndTrim(fns, filts, # fnRs, filtRs,
                     maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=3) # On Windows set multithread=FALSE


# sanity check  comparison of before and after filtration
plotQualityProfile(c(fns[1:2],filts[1:2]))

# LEARN ERROR RATES ####
# Since some samples may have had zero reads pass QC, reassign filts
filts <- sort(list.files(filtpath, full.names = TRUE))
errF <- learnErrors(filts, multithread=TRUE, MAX_CONSIST = 20)

# sanity check for error model
plotErrors(errF, nominalQ=TRUE)

# DEREPLICATION ####
derep <- derepFastq(filts, verbose=TRUE)


# Name the derep-class objects by the sample names
# If some samples were removed (no reads passed QC), reassign sample.names
if(length(derep) != length(sample.names)){
  sample.names <- unlist(map(strsplit(basename(filts), "_filt"), 1))
}
names(derep) <- sample.names

# SAMPLE INFERRENCE ####
dadaFs <- dada(derep, err=errF, multithread=TRUE, selfConsist = TRUE, verbose=TRUE, pool = "pseudo")

# MAKE SEQUENCE TABLE ####
seqtab <- makeSequenceTable(dadaFs)

# REMOVE CHIMERAS ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# reassign "out" to remove any missing reads
out = out[as.data.frame(out)$reads.out > 0,]

# TRACK READS THROUGH PIPELINE ####
getN <- function(x) sum(getUniques(x))
track <- cbind(out[,1], sapply(dadaFs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "nonchim")
rownames(track) <- sample.names
track = as.data.frame(track)
track$total.loss.proportion = (track[,1]-track[,2])/track[,1]
head(track)
# write.csv(track, file = "./output/read_counts_at_each_step.csv", row.names = TRUE)


# Save intermediate seqtable object
saveRDS(seqtab.nochim, "./output/seqtab.nochim.RDS")


# IMPORT METADATA ####
meta = read_delim("./data/Metadata.csv",delim = ",")
row.names(meta) <- as.character(meta$SampleID)
meta <- meta[meta$SampleID %in% sample.names,]

identical(meta$SampleID,row.names(seqtab.nochim))

row.names(meta) <- meta$SampleID


# Remove all seqs with fewer than 100 nucleotides ####
keeper_esvs <- nchar(names(as.data.frame(seqtab.nochim))) > 99
seqtab.nochim <- seqtab.nochim[,keeper_esvs]

# ASSIGN TAXONOMY ####
taxa <- assignTaxonomy(seqtab.nochim, "./taxonomy/UNITE_Euk_2020-02-04_non-dev.fasta.gz", multithread=3)

# Save intermediate files
# write.csv(as.data.frame(seqtab.nochim), file = "./output/SeqTable_no-chimera.csv", row.names = TRUE, quote = FALSE)
saveRDS(seqtab.nochim, file = "./output/dada2_seqtable.RDS")
saveRDS(taxa, file = "./output/RDP_Taxonomy_from_dada2.RDS")

# re-load point
# seqtab.nochim <- readRDS("./output/dada2_seqtable.RDS")
# taxa <- readRDS("./output/RDP_Taxonomy_from_dada2.RDS")
# meta <- read_delim("./data/Metadata.csv",delim = ",")
# row.names(meta) <- as.character(meta$SampleID)


# Hand off to Phyloseq ####
otu <- otu_table(seqtab.nochim,taxa_are_rows = FALSE)
tax <- tax_table(taxa)
met <- sample_data(meta)
row.names(met) <- row.names(meta)


ps <- phyloseq(otu,met,tax)

# Find non-fungi
ps_nonfungi <- subset_taxa(ps, Kingdom != "k__Fungi")

ps %>% transform_sample_counts(function(x){x/sum(x)}) %>%
  plot_bar2(fill="Kingdom")
# ggsave("./output/figs/Kingdom_Level_Taxonomic_Proportions.png",dpi=300)

ps_nonfungi %>% 
  subset_taxa(Kingdom %in% c("k__Metazoa","k__Alveolata","k__Viridiplantae","k__Stramenopila")) %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  plot_bar2(fill="Kingdom")
# ggsave("output/figs/Kingdom_Level_Taxonomic_Proportions_not-including-fungi.png")

ps %>% subset_taxa(Kingdom == "k__Metazoa") %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  plot_bar2(fill="Phylum")
# ggsave("output/figs/Phylum_Level_Taxonomic_Proportions_of_Metazoans.png")

# REMOVE NON-FUNGI and empty samples/taxa ####
ps <- subset_taxa(ps, Kingdom == "k__Fungi")
ps <- subset_taxa(ps, taxa_sums(ps) > 0)
ps <- subset_samples(ps, sample_sums(ps) > 0)

# Save DNA sequences apart from rownames (from subsetted ps object)
seqs <- taxa_names(ps)
seqs <- DNAStringSet(seqs)
saveRDS(seqs,"./output/ASV_reference_sequences.RDS")


pretty_names <- paste("FungalASV",1:length(taxa_names(ps)),":",
      tax_table(ps)[,2],
      tax_table(ps)[,3],
      tax_table(ps)[,4],
      tax_table(ps)[,5],
      tax_table(ps)[,6],
      tax_table(ps)[,7], sep="_") %>%
  str_remove("k__") %>% str_remove("p__") %>% str_remove("c__") %>% str_remove("o__") %>% str_remove("f__") %>% str_remove("g__") %>% str_remove("s__") %>%
  str_replace(pattern = "_:_",replacement = ": ")

df <- data.frame(TaxaName=pretty_names,Sequence=taxa_names(ps))
saveRDS(df,"./output/SequenceNames_and_Taxonomy.RDS")

# Set Seawater as first level of Sponge_Species
ps@sam_data$Sponge_Species <- factor(ps@sam_data$Sponge_Species, 
                                     levels = c("Seawater","Chondrilla","Chondrosia","Crambe","Petrosia"))


# Save RDS object for Phyloseq
saveRDS(ps, file = "./output/clean_phyloseq_object.RDS")


