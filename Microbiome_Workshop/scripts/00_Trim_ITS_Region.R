# -----------------------------------------------------------------------------#
# Microbiome analysis workshop
# Trimming conserved flanking regions away from ITS region
# Author: Geoffrey Zahn
# Software versions:  ITSxpress v 1.8.0
#                     R v 3.6.3
# -----------------------------------------------------------------------------#

#############################################################
#### This script calls itsxpress to remove flanking      ####
#### conserved regions of the rDNA when using fungal     ####
#### ITS1 or ITS2 metabarcodes.                          ####
#### You must have itsxpress installed on your system    ####
#### and present in your PATH. See itsxpress             ####
#### installation documents for instructions.            ####
#############################################################

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")

#################################
#########  PARAMETERS  ##########

# Set the number of parallel threads to use (depends on your hardware):
nthreads <- 4

# Set the region (must be either "ITS1" or "ITS2")
itsregion <- "ITS1"

#################################
#################################

# EXTRACT ITS REGION FROM ALL "cutadapted" SEQUENCES ####

# Run ITSxpress on all forward fastq files (reverse reads often aren't used in fungal analyses)

# find the "cutadapted" files
fwds <- list.files("./data/filtN",pattern = "_pass_1",full.names = TRUE)

# build names for outfiles
outs <- paste0(str_remove_all(fwds,"pass_1.fastq.gz"), "ITS.fastq.gz")

# here's a list of the filenames that will be output:
outs

# build the ITSxpress command and run it on each file in turn
for(i in 1:length(fwds)){
  
  itsxpress <- paste0("itsxpress --fastq ",fwds[i],
                  " --outfile ",outs[i],
                  ".ITS1.gz --region ",itsregion,
                  " --taxa Fungi",
                  " --threads ",nthreads,
                  " --log ",outs[i],".log",
                  " --single-end")

  system(command = itsxpress)
}

# This will generate a new set of output files, tagged with "ITS" in the filenames that can be used for the subsequent 
# filter-and-trim step in our workflow

# an example of the structure of a call to itsxpress is shown below:
"itsxpress --fastq ./data/filtN/2554_pass_1.fastq.gz --outfile ./data/filtN/2554_ITS.fastq.gz.ITS1.gz --region ITS1 --taxa Fungi --threads4 --log ./data/filtN/2554_ITS.fastq.gz.log --single-end"
