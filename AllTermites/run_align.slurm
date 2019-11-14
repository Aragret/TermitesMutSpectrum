#!/bin/bash

#SBATCH --job-name=mafft_mito
#SBATCH --partition=compute
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1
#SBATCH --mail-user=alina.mikhailova@oist.jp
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --cpus-per-task=3
#SBATCH --array=0-12%3
#SBATCH --error=mafft_mito.err

module purge
module load mafft/7.305
module load emboss/6.6.0

cd forAlign/

files=($(ls))

input=${files[${SLURM_ARRAY_TASK_ID}]}

#For protein-coding genes
#Use transeq to get the protein sequences, or use MEGA
transeq  -table 5 -sequence ${input} -outseq results/${input}.pep

#Align protein sequences with MAFFT
mafft --maxiterate 1000 --globalpair results/${input}.pep > results/${input}.align

#Convert your protein alignment in DNA alignment
perl ../pal2nal.pl results/${input}.align ${input} -codontable 5 -output fasta > results/${input}.align.fna
