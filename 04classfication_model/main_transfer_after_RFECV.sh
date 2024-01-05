#!/bin/bash
#SBATCH --job-name=RandomForest+RFECV
#SBATCH --partition=amd
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8000
#SBATCH --nodes=1
#SBATCH --time=90-00:00:00

/home1/jialh/anaconda3/bin/python \
/home1/jialh/brain/01meta/multikingdom/scripts/main_transfer_after_RFECV.py