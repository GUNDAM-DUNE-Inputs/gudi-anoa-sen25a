#!/usr/bin/env bash
#SBATCH -p extended-96core-shared
#SBATCH --array=0-42                           # 6 params × (2·NUM_SIGMA+1)=5 (NUM_SIGMA = 2) σ-points ⇒ 30 tasks
#SBATCH --cpus-per-task=8
#SBATCH -t 01:00:00
#SBATCH --mem=50G
set -euo pipefail
set -x

# ── Config ────────────────────────────────────────────────────────────────────
input_yaml="./overrides/MaCh3Comparison300Zx300Esmoothing04Profiling.yaml"

# parameter definitions (must stay in sync)
paramName=(sin12 sin13 sin23 m21 m32 dcp)
occurrence=(1    2     3     4   5    6  )
priorDefault=(0.303 0.02225 0.452 7.41e-5 2.51e-3 -2.233)
sigma=(        0.013  0.0007   0.021 1.8e-6     3.4e-5     1.17   )
#sigma=(        0.013  0.0007   0.021 1.8e-6     1.05e-4     0.316   ) #- 10 points
#sigma=(        0.013  0.0007   0.0217 1.8e-6     1.41e-4     1.09 ) #- 7 points

NUM_SIGMA=3                                    # ±2σ ⇒ 5 points (−2,−1,0,+1,+2)
STEP_COUNT=$(( 2*NUM_SIGMA + 1 ))               # =5
PARAM_COUNT=${#paramName[@]}                   # =6

# ── Map SLURM task → (param index, σ-offset) ─────────────────────────────────
TASK_ID=$SLURM_ARRAY_TASK_ID                    # 0…29
param_idx=$(( TASK_ID / STEP_COUNT ))           # integer division
step_idx=$(( TASK_ID % STEP_COUNT ))            # 0…4
offset=$(( step_idx - NUM_SIGMA ))              # maps 0→−2,1→−1,2→0,3→+1,4→+2

# ── Pull out the right values ────────────────────────────────────────────────
name=${paramName[param_idx]}
occ=${occurrence[param_idx]}
default=${priorDefault[param_idx]}
step=${sigma[param_idx]}

# ── Compute the new prior = default + offset·σ ───────────────────────────────
newPrior=$(awk -v d=$default -v s=$step -v k=$offset \
  'BEGIN{ printf "%g\n", d + k*s }')

echo "occ=$occ, prior=$newPrior, name=$name"

echo "Task#$TASK_ID → $name (occ=$occ), offset=$offset σ ⇒ prior=$newPrior"
#echo "./submit_DUNE_profiling.sh "$occ" "$newPrior" "$name" "
# ── Launch your profiling script ─────────────────────────────────────────────
./submit_DUNE_profiling.sh "$occ" "$newPrior" "$name" "$input_yaml"