#!/bin/bash

#SBATCH --job-name=iqtree_mito
#SBATCH --partition=compute
#SBATCH --time=50:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1
#SBATCH --mail-user=alina.mikhailova@oist.jp
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --cpus-per-task=12
#SBATCH --error=iqtree_mito.err


./src/iqtree -s align_cat.fasta -nt 12 -m TEST -m MFP+MERGE -rcluster-max 2000 -rcluster 10 -alrt 1000 -bb 1000 -bnni -pre ${SLURM_JOB_NAME} -asr

