#!/bin/bash
#SBATCH --time=300 # days-hours:minutes:seconds (or just minutes) ... how much time are you asking for (max 3 days)
#SBATCH --nodes=1 # How many nodes are you asking for? Unless you know a lot about parallel computing, stick to 1     
#SBATCH -o slurm-%j.out # the name of your STDOUT (this adds the jobid and can stay the same for every job)
#SBATCH -e slurm-%j.err # the name of your STDERR (this adds the jobid and can stay the same for every job)
#SBATCH --ntasks=16 # how many cores are you asking for?
#SBATCH --account=uvu # what is your account (students will put "uvu-biol490r" here)
#SBATCH --partition=kingspeak-shared # which partition do you want?
#SBATCH --mail-type=END,FAIL # should they email you? when?               
#SBATCH --mail-user=your_email@uvu.edu # what email address do you want updates sent to? Change this, obviously.
#SBATCH --mem=64000 # How much memory are you requesting (in Mb)

# load any modules your script needs
module load blast

# set up any shell variables you need
QUERY=/scratch/general/vast/Zahn/Mangrove/Mangrove.fasta
BLAST_OUT=~/Mangrove_blast_out.txt

# your actual job:

blastn -query $QUERY -db NCBI_VIRUS -outfmt '6 qseqid sseqid sskingdom ssciname sblastname pident evalue length bitscore gaps' -out $BLAST_OUT -num_threads 16 -max_target_seqs 1