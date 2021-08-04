# -----------------------------------------------------------------------------#
# Microbiome analysis workshop
# Processing raw amplicon reads
# Author: Geoffrey Zahn
# Software versions:  R v 3.6.3
#                     tidyverse v 1.3.0
#                     dada2 v 1.14.1
#                     decontam v 1.6.0
#                     phyloseq v 1.30.0
#                     purrr v 0.3.4
#                     Biostrings v 2.54.0
#                     patchwork v 1.0.1
# -----------------------------------------------------------------------------#

#################################################################################
#                               Main workflow                                   #
# Filter and trim, denoise, sample inferrence, chimera and contaminant removal, # 
# taxonomic assignment, combine sequence table and metadata                     #
#################################################################################

# PACKAGES, SCRIPTS, AND SETUP ####

# why each package (put in onboarding document)
library(tidyverse); packageVersion("tidyverse")
library(dada2); packageVersion("dada2")
library(decontam); packageVersion("decontam")
library(phyloseq); packageVersion("phyloseq")
library(purrr); packageVersion("purrr")
library(Biostrings); packageVersion("Biostrings")
library(patchwork); packageVersion("patchwork")



# PARSE FILE PATHS ####

# File parsing - 

path <- "./data/filtN" # CHANGE to the directory containing your demultiplexed fastq files when using your own data
filtpath <- file.path(path, "filtered") # Filtered files go into the filtered/ subdirectory
if(!file_test("-d", filtpath)) dir.create(filtpath) # make directory for filtered fqs if not already present


# Here, the pattern is set to match the forward reads tagged as "...pass_1.fastq.gz" in their filenames; 
# "...pass_2.fastq.gz" for reverse reads

# Your data may differ, using "F" and "R" in the filenames, or something similar..
# Be sure to change that pattern to fit your own files when using your own data
fns <- sort(list.files(file.path(path), full.names = TRUE, pattern = "pass_1.fastq.gz")) # make pattern match your FWD reads
rns <- sort(list.files(file.path(path), full.names = TRUE, pattern = "pass_2.fastq.gz")) # make pattern match your REV reads

# This may have to be modified for your own data if your files follow a differnet naming convention
# The following line splits the filename on "_" and keeps the first element
sample.names <- unlist(map(strsplit(basename(fns), "_"), 1)) # this pulls out just the basename

# visualize a couple of fwd read quality profiles to help select reasonable filtration parameters
# you can select any number of files here...
# as-is, this just shows the Fwd and Rev quality profiles for the 1st and 2nd files
p1 <- plotQualityProfile(fns[1:2]) + ggtitle("Example forward reads")
p2 <- plotQualityProfile(rns[1:2]) + ggtitle("Example reverse reads")

# display and save the plots
p1 / p2
ggsave("./output/figs/unfiltered_quality_plots.png",dpi=500,height = 6,width = 6)

# FILTER AND TRIM ####

# here, we decide what filenames we want our filtered reads to be called
# in this case, we use the base name of the sample and save the filtered files into the "filtered" subdirectory
filts_f <- file.path(path, "filtered", paste0(sample.names, "_FWD_filt.fastq.gz"))
filts_r <- file.path(path, "filtered", paste0(sample.names, "_REV_filt.fastq.gz"))

# this is the actual qualit control step
# These values are informed by our quality plots
out <- filterAndTrim(fns, filts_f, rns, filts_r, # input and output file names as denoted above
                     maxN=0, # uncalled bases are currently not supported in dada2
                     maxEE=c(2,2), # refers to the maximum expected errors allowed
                     truncQ=2, # special value denoting "end of good quality sequence" (optional)
                     rm.phix=TRUE, # automatically remove PhiX spike-in reads from sequencing center
                     truncLen = c(250,200), # refers to the lengths at which to truncate Fwd and Rev reads, respectively 
                     compress=TRUE, # compress output files with gzip
                     multithread=4) # On Windows set multithread=FALSE

# save the filtration info in case we need it later
saveRDS(out, "./output/trackreads.RDS")

# Did any samples have NO reads pass the filtration step?

length(fns);length(filts_f) # will be the same length if all samples had some passing reads
length(rns);length(filts_r) # will be the same length if all samples had some passing reads

# In case some samples may have had zero reads pass QC, reassign filts
filts_f <- sort(list.files(filtpath, full.names = TRUE,pattern = "FWD"))
filts_r <- sort(list.files(filtpath, full.names = TRUE,pattern = "REV"))

# sanity check  comparison of before and after filtration
p3 <- plotQualityProfile(fns[1]) + ggtitle("Unfiltered")
p4 <- plotQualityProfile(filts_f[1])+ ggtitle("Filtered") + coord_cartesian(xlim = c(0,300))
p3 / p4
ggsave("./output/figs/filtered_quality_comparison.png",dpi=500,height = 6,width = 6)

# LEARN ERROR RATES ####

# learn errors
set.seed(123) # "random" seed for reproducibility
errF <- learnErrors(filts_f, multithread=TRUE, MAX_CONSIST = 20,verbose = 2,randomize = TRUE) # set multithread = FALSE on Windows
saveRDS(errF,"./output/errF.RDS")
set.seed(123)
errR <- learnErrors(filts_r, multithread=TRUE, MAX_CONSIST = 20,verbose = 2,randomize = TRUE) # set multithread = FALSE on Windows
saveRDS(errR,"./output/errR.RDS")

# sanity check for error model
# explain what to look for in the plot! 
plotErrors(errF, nominalQ=FALSE)
ggsave("./output/figs/error_model.png",dpi=500,height = 6,width = 6)
plotErrors(errR, nominalQ=FALSE)

# dereplication
derepF <- derepFastq(filts_f, verbose=TRUE)
derepR <- derepFastq(filts_r, verbose=TRUE)

# Name the derep-class objects by the sample names
# If some samples were removed (no reads passed QC), reassign sample.names
if(length(derepF) != length(sample.names)){
  sample.names <- unlist(map(strsplit(basename(filts_f), "_filt"), 1))
}


if(identical(unlist(map(strsplit(basename(filts_f), "FWD_filt"), 1)),unlist(map(strsplit(basename(filts_r), "REV_filt"), 1)))){
  names(derepF) <- sample.names
  names(derepR) <- sample.names
} else {
  stop("Make sure fwd and rev files are in same order!")
}  



# SAMPLE INFERRENCE ####
dadaFs <- dada(derepF, err=errF, multithread=TRUE, selfConsist = TRUE, verbose=TRUE, pool = "pseudo") # set multithread = FALSE on Windows
dadaRs <- dada(derepR, err=errR, multithread=TRUE, selfConsist = TRUE, verbose=TRUE, pool = "pseudo") # set multithread = FALSE on Windows

# MERGE FWD and REV READS ####
mergers <- mergePairs(dadaFs, filts_f, dadaRs, filts_r, verbose=TRUE)

# MAKE SEQUENCE TABLE ####
seqtab <- makeSequenceTable(mergers)

# REMOVE CHIMERAS ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
saveRDS(seqtab.nochim,"./output/seqtab.nochim.RDS")
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
?removeBimeraDenovo
chimera.time <- Sys.time()

# reassign "out" to remove any missing reads
out = out[as.data.frame(out)$reads.out > 0,]

# TRACK READS THROUGH PIPELINE ####
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track = as.data.frame(track)
track$total.loss.proportion = (track[,1]-track$nonchim)/track[,1]
head(track)

write.csv(track, file = "./output/16S_read_counts_at_each_step.csv", row.names = TRUE)


readLines("./output/16S_read_counts_at_each_step.csv")

# IMPORT METADATA ####
meta <- read_csv("./metadata/chagos_metadata_clean.csv")
row.names(meta) <- as.character(meta$Sample_ID)

# reorder metadata to match seqtab
df <- data.frame(seqtab_rows=row.names(seqtab.nochim),
                 Sample_ID=row.names(seqtab.nochim))
df2 <- left_join(meta,df,by="Sample_ID")
row.names(df2) <- df2$SRA_Accession
meta <- df2[row.names(seqtab.nochim),]
row.names(meta) <- meta$SRA_Accession
identical(row.names(meta),row.names(seqtab.nochim))


# Remove all seqs with fewer than 100 nucleotides (if any) ####
keeper_esvs <- nchar(names(as.data.frame(seqtab.nochim))) > 99
seqtab.nochim <- seqtab.nochim[,keeper_esvs]

# remove singleton taxa
seqtab.nochim <- seqtab.nochim[,colSums(seqtab.nochim) > 1]

# remove newly singleton taxa
seqtab.nochim <- seqtab.nochim[,colSums(seqtab.nochim) > 1]

# save cleaned up seqtab
saveRDS(seqtab.nochim,"./output/seqtab.nochim.clean.RDS")

# Find and remove contaminants ####

# need to explain how to identify these in data set
meta$Control <- meta$Host_ID == "Blank"
meta <- meta[meta$SRA_Accession %in% row.names(seqtab.nochim),] 

contams = isContaminant(seqtab.nochim, neg = meta$Control, normalize = TRUE)
library(decontam)
?isContaminant
table(contams$contaminant) # how many taxa are contaminants?
write.csv(contams, file = "./output/likely_contaminants.csv", row.names = TRUE)

# remove contaminant sequences and control samples from both tables, respectively ####
seqtab.nochim = seqtab.nochim[,(which(contams$contaminant != TRUE))]
seqtab.nochim = seqtab.nochim[meta$Control == FALSE,]
meta = meta[meta$Control == FALSE,]


# ASSIGN TAXONOMY ####

# need to explain how this is done (databases, etc.)
taxa <- assignTaxonomy(seqtab.nochim, "./taxonomy/rdp_train_set_16.fa.gz", multithread=4)

# Save intermediate taxonomy file
saveRDS(taxa, file = "./output/RDP_Taxonomy_from_dada2.RDS")

# add_species names
taxa <- addSpecies(taxa, "./taxonomy/rdp_species_assignment_16.fa.gz")

# Save completed taxonomy file
saveRDS(taxa, file = "./output/RDP_Taxonomy_from_dada2_sp.RDS")


# inspect taxonomy
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


# Hand off to Phyloseq ####
otu <- otu_table(seqtab.nochim,taxa_are_rows = FALSE)
tax <- tax_table(taxa)
met <- sample_data(meta)
row.names(met) <- meta$SRA_Accession


ps <- phyloseq(otu,met,tax)
saveRDS(ps,"./output/ps_not-cleaned.RDS")




