#!/bin/bash
#PBS -j oe
#PBS -N blastAMPstoproteome
#PBS -l walltime=2:00:00
#PBS -l select=1:ncpus=28:mem=32gb
#PBS -m ae

cd $PBS_O_WORKDIR
shopt -s expand_aliases

module load conda3
source $CONDA_PROF/conda.sh
conda activate blast-2.11.0

gzip -d M_musculus-proteome-UP000000589.fasta.gz 

makeblastdb -in Mus_musculus.fasta -dbtype 'prot'

blastp -db Mus_musculus.fasta -query M_musculus-proteome-UP000000589.fasta -outfmt 6 -max_target_seqs 5 -evalue=10 -num_threads 10 > Mus_musculus_proteome.blastp


rm M_musculus-proteome-UP000000589.fasta
rm Mus_musculus.fasta.* 