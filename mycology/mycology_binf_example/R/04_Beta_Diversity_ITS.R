# -----------------------------------------------------------------------------#
# Mediterranean sponge fungi - ocean acidification project 
# Beta-Diversity measures nd Ordinations
# Author: Geoffrey Zahn
# Software versions:  R v 3.6.3
#                     tidyverse v 1.3.0
#                     patchwork v 1.0.1
#                     phyloseq v 1.30.0
#                     broom v 0.7.0
#                     purrr v 0.3.4
#                     vegan v 2.5.6
#                     ade4 v 1.7.15
#                     microbiome v 1.8.0
# -----------------------------------------------------------------------------#

# packages ####
library(tidyverse); packageVersion("tidyverse")
library(patchwork); packageVersion("patchwork")
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(broom); packageVersion("broom")
library(purrr); packageVersion("purrr")
library(ade4); packageVersion("ade4")
library(vegan); packageVersion("vegan")
# BiocManager::install("microbiome")
library(microbiome); packageVersion("microbiome")


#functions
source("./R/palettes.R")


# IMPORT DATA ####
ps <- readRDS("./output/clean_phyloseq_object.RDS")

# remove "NA" Phylum taxa
ps <- subset_taxa(ps,!is.na(tax_table(ps)[,2]))

# remove empty sites
ps <- ps %>% subset_samples(sample_sums(ps) > 0)

# relative abundance transformation
ps_ra <- ps %>% transform_sample_counts(function(x){x/sum(x)}) %>% subset_samples(Sponge_Species != "Seawater")

# Ordinations ####
set.seed(123)
ord <- ordinate(ps_ra,method = "DCA",distance = "jaccard")
p1 <- plot_ordination(ps_ra,ord,color="Sampling_Site") + scale_color_manual(values=pal.discrete) + theme(legend.position = "top") +
  labs(color="Site")
p2 <- plot_ordination(ps_ra,ord,color="Sponge_Species") + scale_color_manual(values=pal.discrete) + 
  theme(legend.position = "top",
        legend.text = element_text(face="italic")) + labs(color="Sponge species")
p3 <- plot_ordination(ps_ra,ord,color="Acidified") + scale_color_manual(values=pal.discrete) + theme(legend.position = "top") +
  labs(color="Acidification status")

p1+p2+p3


# Beta-dispersion
names(meta(ps_ra))
w <- betadiver(otu_table(ps_ra),"w")
w.disper <- betadisper(w,group = meta(ps_ra)$Sponge_Species)

plot(w.disper)
