#!/bin/bash

#PBS -q eternity ## big eternity mem1t mem512g
#PBS -d .
#PBS -l walltime=300:00:00,mem=10G,nodes=1:ppn=4

autoPipeline_simple.py COIIRefSeqs.fasta Cox2 > /mnt/lustre/agmikhaylova/termits/2Derived/result.out


