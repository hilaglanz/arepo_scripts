#!/bin/bash

#SBATCH --time=3-0
#SBATCH -J arepo_plotting

### Output files 
#SBATCH --output=/storage/ph_hagai/glanz/ThickAcretion/plots/job.%J.out.txt
#SBATCH --error=/storage/ph_hagai/glanz/ThickAcretion/plots/job.%J.err.txt

export SOURCE_DIR = "/storage/ph_hagai/glanz/ThickAcretion/output"
export SAVING_DIR = "/storage/ph_hagai/glanz/ThickAcretion/plots"

arepoPy make_plots.py --source_dir=$SOURCE_DIR --saving_dir=$SAVING_DIR



