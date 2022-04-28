#!/bin/bash

#SBATCH -n 1
#SBATCH --time=3-0
#SBATCH -J arepo_plotting

### Output files 
#SBATCH --output=/storage/ph_hagai/glanz/ThickAcretion/plots/job.%J.out.txt
#SBATCH --error=/storage/ph_hagai/glanz/ThickAcretion/plots/job.%J.err.txt

export SOURCE_DIR = "/storage/ph_hagai/glanz/ThickAcretion/output"
export SAVING_DIR = "/storage/ph_hagai/glanz/ThickAcretion/plots"

export arepoPy = "/usr/local/ph_hagai/anaconda3/envs/amuse-env/bin/python-3.10"

$arepoPy make_plots.py --source_dir=$SOURCE_DIR --saving_dir=$SAVING_DIR



