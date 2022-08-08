#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate snakemake 

snakemake -c 72 --use-singularity --singularity-args "-B /projects,/home,/gsc" --rerun-incomplete --restart-times=3
