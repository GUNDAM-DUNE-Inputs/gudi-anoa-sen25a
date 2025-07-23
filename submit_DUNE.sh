#!/bin/bash
#SBATCH -p extended-96core-shared
#SBATCH --output=job_output_%A_%a.out                       
#SBATCH --cpus-per-task=16
#SBATCH -t 24:00:00 
#SBATCH --mem=35G

## Set Up the GUNDAM Environment
# at seawulf
source /gpfs/home/uyevarouskay/Work/env_gundam_home_v4.sh

# at SBU NN home cluster
#source /home/isgould/work/gundam/install/setup.sh

## Set the Input Files Path
# at seawulf cluster
export OA_INPUT_FOLDER=/gpfs/scratch/uyevarouskay/atm/gundam_files_fin_v13/ 

# at SBU NN home cluster
#export OA_INPUT_FOLDER=/storage/shared/DUNE/OA-inputs/atm/gudi-inputs/v3/

# dunegpvm
#export OA_INPUT_FOLDER=/pnfs/dune/persistent/users/weishi/OA-inputs/atm_reprocessed_v2/ #might be outdated

## NuOscillator 
export NUOSCILLATOR_ROOT_LIB=./gundamOscAnaTools/resources/TabulateNuOscillator/build-x86_64/
source ${NUOSCILLATOR_ROOT_LIB}/bin/setup.NuOscillator.sh

# if you request an Asimov fit, simply disregard the `-d` option, --scan is to build llh scans
gundamFitter -a -d --scan -c ./config_DUNE.yaml -t 8
# To disable  oscillation parameters consider override the parameters config:
#gundamConfigUnfolder -c config_DUNE.yaml -of overrides/disableOscillationParameters.yaml -o
#and run instead
#gundamFitter -a -d -c ./config_DUNE_With_disableOscillationParameters.json -t 8
