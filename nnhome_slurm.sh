#!/bin/bash
#SBATCH --job-name=ATM
#SBATCH --output=slurm_files/GundamInputDUNE_LBL_%j.out
#SBATCH --error=slurm_files/GundamInputDUNE_LBL_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=[youremailhere@yourschool.edu]
#SBATCH --mem=30G
#SBATCH -c 16
#SBATCH --gres=gpu

#source /home/rrazakami/workspace/ROOT/root_binary/bin/thisroot.sh
source /home/isgould/work/gundam/install/setup.sh
#source /home/fyguo/GUNDAM/Install/gundam/setup.sh

#display detailed information about your system
#uname -a
#display the date
#date
#echo "started"

#run fitter
#echo "gundamFitter -c config_DUNE.yaml -d -t 16"
#gundamFitter -c config_DUNE.yaml -d -t 16
# gundamFitter -c config_DUNE.yaml -a --scan -t 16
# time gundamFitter -a -c config_DUNE.yaml --gpu -t 16
# time gundamFitter -a -c config_DUNE.yaml --gpu --kick-mc 0.1 -t 16
# time gundamFitter -a -c config_DUNE.yaml --gpu --kick-mc 1 -t 16
# time gundamFitter -a -c config_DUNE.yaml --gpu -d --toy -t 16

./submit_DUNE.sh

#echo "finished"
#date
