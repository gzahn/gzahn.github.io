



# PREPARE YOUR SYSTEM TO INSTALL STUFF ####

# Highlight everything in this script (CTL-A) and then run it (CTL-ENTER)



# install devtools package which allows you access to other tools
install.packages("devtools")

# install Bioconductor (set of bioinformatics tools)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.13")

# "STANDARD" INSTALLATIONS ####

install.packages("ape", dependencies = TRUE )
install.packages("corncob", dependencies = TRUE )
install.packages("msa", dependencies = TRUE )
install.packages("patchwork", dependencies = TRUE )
install.packages("phangorn", dependencies = TRUE )
install.packages("seqinr", dependencies = TRUE )
install.packages("tidyverse", dependencies = TRUE )
install.packages("vegan", dependencies = TRUE )

# "NON-STANDARD" INSTALLATIONS ####

# use bioconductor to install some specific and useful tools
BiocManager::install("Biostrings")
BiocManager::install("ShortRead") 
BiocManager::install("decontam")
BiocManager::install("phyloseq")

# use devtools to install dada2 package
devtools::install_github("benjjneb/dada2", ref="v1.18")


# test to make sure package installation worked ####
library(ape)
library(corncob)
library(msa)
library(patchwork)
library(phangorn)
library(tidyverse)
library(vegan)
library(Biostrings)
library(ShortRead)
library(decontam)
library(phyloseq)
library(dada2)

packages <- c("ape","corncob","msa","patchwork","phangorn","tidyverse","vegan",
              "Biostrings","ShortRead","decontam","phyloseq","dada2")
'%ni%' <- Negate("%in%")

not_loaded <- packages[(packages %ni% (.packages()))]

if(length(not_loaded) > 0){
  paste0("The following package is not loaded properly: ",not_loaded)
} else {
  print("All the packages were loaded properly")
}

