#!/bin/bash
#SBATCH -p extended-96core-shared
#SBATCH --output=job_output_%A_%a.out                       
#SBATCH --cpus-per-task=1
#SBATCH -t 48:00:00 
#SBATCH --mem=10G

module load root

python example.py
