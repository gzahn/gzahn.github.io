#!/bin/bash
set -ueo pipefail

Rscript ./scripts/01_Process_Raw_16S_Reads.R
Rscript ./scripts/02_Build_and_add_Phylogeny.R
Rscript ./scripts/03_Clean_and_Explore_Data.R
Rscript ./scripts/04_Explore_and_Test_Hypotheses.R

