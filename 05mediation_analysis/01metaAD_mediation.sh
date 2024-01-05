#!/bin/bash
#SBATCH --job-name=mediation
#SBATCH --partition=DCU
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8000
#SBATCH --nodes=1
#SBATCH --time=9000-00:00:00

#/home1/jialh/anaconda3/envs/R420/bin/Rscript \
/home1/jialh/tools/anaconda3/bin/python \
/home1/jialh/brain/01meta/multikingdom/04mediation/01metaAD_mediationAnalysis_parallel.py