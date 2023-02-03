#!/bin/bash

# Get all argument and put it in correct format
sampleName=$(echo $@ | tr -s '[:blank:]' ',')
sampleName="[${sampleName}]"

eval "$(conda shell.bash hook)"
conda activate snakemake 

snakemake -c 72 --use-singularity --singularity-args "-B /projects,/home,/gsc" --rerun-incomplete --restart-times=3 --config sample=$sampleName 
