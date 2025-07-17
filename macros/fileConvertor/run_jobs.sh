#!/bin/bash

# On SeaWulf
 module load root
# On NNhome machine
#source /home/rrazakami/workspace/ROOT/root_binary/bin/thisroot.sh

# Adjust for the local data location.
export OA_INPUTS=/gpfs/scratch/uyevarouskay/atm/MaCH3-inputs/
export OA_OUTPUTS=/gpfs/scratch/uyevarouskay/atm/gundam_files_fin_v11/

# Count the number of .root files
num_files=$(ls ${OA_INPUTS}/*.root | wc -l)

# Loop over numeric indices: 0 to num_files-1
for ((i=0; i<num_files; i++)); do
    echo "./fileConvertor_ATM $i"
    ./fileConvertor_ATM "$i"
done