#!/bin/bash

#SBATCH --job-name=Fin16b
#SBATCH --partition=compute
#SBATCH --time=1-0
#SBATCH --mem-per-cpu=5G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mail-user=alina.mikhailova@oist.jp
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --input=none
#SBATCH --output=job_%j.out

InputAlignment=FcC_supermatrix.fas_with3rd
PartitioningScheme=geneAndCodon_partitinoning_scheme_with3rd.txt

mkdir ${SLURM_JOB_NAME}
cd ${SLURM_JOB_NAME}

../../src/iqtree -s ../${InputAlignment} -nt ${SLURM_CPUS_PER_TASK} -spp ../${PartitioningScheme} -mset GTR -mrate I,G,I+G -m TESTMERGE -rcluster-max 2000 -rcluster 10 -alrt 1000 -bb 1000 -bnni -pre ${SLURM_JOB_NAME} -asr -o Empusa_pennata,Mantis_religiosa,Metallyticus_splendidus
