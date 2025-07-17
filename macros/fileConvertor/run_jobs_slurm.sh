#!/bin/bash
#SBATCH -p extended-40core-shared
#SBATCH --job-name=gundam-inputs-atm
#SBATCH --output=slurm_files/slurm-%A_%a.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --array=1-36          # One task per file
#SBATCH --time=02:00:00

# Load ROOT module
module load root

# Set paths
export OA_INPUTS=/gpfs/scratch/uyevarouskay/atm/MaCH3-inputs/
export OA_OUTPUTS=/gpfs/scratch/uyevarouskay/atm/gundam_files_fin_v12/

# Read all files into an array
files=()
for file in ${OA_INPUTS}/*; do
    files+=("$file")
done

# Calculate correct index (0-based for bash arrays)
index=$((SLURM_ARRAY_TASK_ID - 1))
selected_file=${files[$index]}

# Run the converter with the correct index (pass index to your C++ program)
echo "Running ./fileConvertor_ATM $index for file: $selected_file"
./fileConvertor_ATM "$index"
