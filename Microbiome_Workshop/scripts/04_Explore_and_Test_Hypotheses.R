# -----------------------------------------------------------------------------#
# Microbiome analysis workshop
# Exploring cleaned data using phyloseq and corncob packages
# Author: Geoffrey Zahn
# Software versions:  R v 3.6.3
#                     tidyverse v 1.3.0
#                     phyloseq v 1.30.0
#                     purrr v 0.3.4
#                     Biostrings v 2.54.0
#                     corncob v 0.1.0
#                     vegan v 2.6.0
#                     patchwork v 1.0.1
# -----------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(purrr); packageVersion("purrr")
library(Biostrings); packageVersion("Biostrings")
library(corncob); packageVersion("corncob")
library(vegan); packageVersion("vegan")
library(patchwork); packageVersion("patchwork")

source("./scripts/bbdml_helper.R")

#################################################################################
#                               Main workflow                                   #
#  Explore alpha and beta diversity, visualize data set, test hypotheses,       #
#  search for differentially abundant taxa                                      #
#                                                                               #
#################################################################################


# Load cleaned phyloseq object ####

ps <- readRDS("./output/clean_phyloseq_object.RDS")


# Getting to know your phyloseq data ####

# number of taxa
ntaxa(ps)

# number of samples
nsamples(ps)

# sample names
sample_names(ps)
rank_names(ps)
# taxa names
taxa_names(ps)

# overview of taxonomy
tax_table(ps)
tax_table(ps)[,1] %>% table()
tax_table(ps)[,2] %>% table()
tax_table(ps)[,3] %>% table()
tax_table(ps)[,4] %>% table()
tax_table(ps)[,5] %>% table()
tax_table(ps)[,6] %>% table()

# sample metadata
sample_data(ps) %>% as_tibble()

# ASV table
otu_table(ps)

# how many sequences observed in each sample?
otu_table(ps) %>% rowSums()

# how many times was each taxon observed?
otu_table(ps) %>% colSums()

# how many different samples was each taxon found in?
asv <- otu_table(ps) %>% as("matrix") %>% as.data.frame() # convert to matrix before you can convert to data frame
asv[asv>0] <- 1 # convert to presence/absence
colSums(asv) # sum of presence (present = 1; absent = 0)

# what was the most widespread taxon (not abundance)
widespread_seq <- names(asv)[which(colSums(asv) == max(colSums(asv)))] # this gives long sequence
tax_table(ps)[widespread_seq,] # this pull that row from ASV table

# what was most abundant (raw counts) taxon?
abund_seq <- which(otu_table(ps) %>% colSums() == max(otu_table(ps) %>% colSums()))
tax_table(ps)[abund_seq,]
otu_table(ps)[,abund_seq]

# access the phylogenetic tree
phy_tree(ps)
plot_tree(ps,color="Class")
# this tree sucks, but it's fast
# explain where to update tree


# Alpha diversity metrics ####

# get table of alpha diversity measures
estimate_richness(ps) # explain output
names(sample_data(ps))

# since we have no singletons, we cannot use Chao1 estimate!
# plot alpha diversity for every sample
plot_richness(ps, 
              measures = c("Observed","Shannon","Simpson"), 
              color = "Colony_Color", 
              sortby = "Observed") +
  theme(axis.text.x = element_blank())

# plot, grouped by colony color with added boxplot
plot_richness(ps, 
              x = "Colony_Color",
              measures = c("Observed","Shannon","Simpson"), 
              color = "Colony_Color", 
              sortby = "Observed") +
  geom_boxplot(alpha = .5) +
  theme_minimal()
ggsave("./output/figs/alpha_diversity_boxplot.png",dpi=300,height = 4,width = 6)

# transform raw counts to relative abundance ####
ps_ra <- transform_sample_counts(ps, fun = function(x){x/sum(x)})

# Beta-diversity ####

# Ordination
dca <- ordinate(ps_ra)
plot_ordination(ps_ra,dca,color = "Colony_Color") 

(
  ord1 <- plot_ordination(ps_ra,dca,color = "Colony_Color",shape="Island") +
  geom_point(size=4)  + theme_minimal() +
    theme(legend.position = "top") +
    labs(title = "DCA - Bray")
)

# try another ordination method
nmds <- ordinate(ps_ra,method = "NMDS")

(
ord2 <- plot_ordination(ps_ra,nmds,color = "Colony_Color",shape="Island") +
  geom_point(size=4)  + theme_minimal() +
    theme(legend.position = "none") +
    labs(title = "NMDS - Bray")
)


# also try with unifrac distance, which takes phylogeny into account

# pull out components
asv <- otu_table(ps_ra) %>% as("matrix") %>% as.data.frame()
meta <- sample_data(ps_ra) %>% as.data.frame()


unifrac.dist <- UniFrac(ps_ra)
betadisper(unifrac.dist,group = meta$Colony_Color) %>% plot() # plot beta-dispersion


unifrac <- ordinate(ps_ra,method = "NMDS",distance = unifrac.dist)

(
ord3 <- plot_ordination(ps_ra,unifrac,color = "Colony_Color",shape="Island") +
  geom_point(size=4) + theme_minimal() +
    theme(legend.position = "none") +
    labs(title = "NMDS - Unifrac")
)

# combine all plots into one figure for comparison
ord1 / ord2 / ord3
ggsave("./output/figs/ordinations.png",dpi=300,width = 6,height = 8)

# permanova ####

# run permanova model with colony_color and Island as predictors (with interaction term included)
permanova.bray <- vegan::adonis(asv ~ meta$Colony_Color * meta$Island,method = "bray")
permanova.bray %>% broom::tidy() %>% View

# try with jaccard distance as well
permanova.jaccard <- vegan::adonis(asv ~ meta$Colony_Color * meta$Island,method = "jaccard")
permanova.jaccard


# Differential abundance/dispersion tests ####

# use non-transformed data!
set.seed(123)
da_analysis_colcolor <- differentialTest(formula = ~ Colony_Color, #abundance
                                             phi.formula = ~ 1, #dispersion
                                             formula_null = ~ 1, #mean
                                             phi.formula_null = ~ 1,
                                             test = "Wald", boot = FALSE,
                                             data = ps,
                                             fdr_cutoff = 0.05,
                                             full_output = TRUE)
plot(da_analysis_colcolor)


# find the significant taxa
da_analysis_colcolor$significant_taxa
da_analysis_colcolor$significant_taxa %>% otu_to_taxonomy(data=ps)

# This is a helper function I wrote. It's found in "scripts/bbdml_helper.R" 
bbdml_obj <- multi_bbdml(da_analysis_colcolor,
                         ps_object = ps,
                         mu_predictor = "Colony_Color",
                         phi_predictor = "Colony_Color",
                         taxlevels = 6)
ps %>% sample_data()
# another helper function found in the same file
plot_multi_bbdml(bbdml_obj,
                 color="Colony_Color", 
                 pointsize = 3)


# This saves a plot for each significant taxon in your environment called bbdml_plot_N (1 through number of sig. taxa)

# view all the plots together using patchwork package
(bbdml_plot_1 / bbdml_plot_2 / bbdml_plot_3) | (bbdml_plot_4 / bbdml_plot_5 / bbdml_plot_6)
ggsave("./output/figs/bbdml_plot_all_sig_taxa.png",height = 6, width = 12)

