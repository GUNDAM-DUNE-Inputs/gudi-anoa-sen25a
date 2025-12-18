#!/bin/bash

## Set Up the GUNDAM Environment
# at seawulf
source /gpfs/home/uyevarouskay/Work/env_gundam_home_v7.sh; #latest version

# at SBU NN home cluster
#source /home/isgould/work/gundam/install/setup.sh

## Set the Input Files Path
# at seawulf cluster
export OA_INPUT_FOLDER=/gpfs/scratch/uyevarouskay/atm/v3/ 

# at SBU NN home cluster
#export OA_INPUT_FOLDER=/storage/shared/DUNE/OA-inputs/atm/gudi-inputs/v3/

# dunegpvm
#export OA_INPUT_FOLDER=/pnfs/dune/persistent/users/weishi/OA-inputs/atm_reprocessed_v2/ #might be outdated

## NuOscillator 
export NUOSCILLATOR_ROOT_LIB=../../gundamOscAnaTools/resources/TabulateNuOscillator/build-x86_64/
source ${NUOSCILLATOR_ROOT_LIB}/bin/setup.NuOscillator.sh

#!/usr/bin/env bash
set -euo pipefail
set -x

occurrence="${1:?Need which isFrozen occurrence to flip (1..6))}"
newPrior="${2:?Need new priorValue (e.g. 0.403)}"
paramName="${3:?Need parameter basename (e.g. sin23_0)}"
input_yaml="${4:?Need a parameter config (e.g. ./overrides/MaCh3Comparison300Zx300Esmoothing04Profiling.yaml)}"

#paramName="PMNS_SIN_SQUARED_23"
#occurrence=3 #(3rd in the list of the parameters)
#newPrior=0.403
#input_yaml="./overrides/MaCh3Comparison300Zx300Esmoothing04Profiling.yaml"
#input_yaml="./overrides/MaCh3ComparisonFixedH300Zx300Eprofiling.yaml"
#input_yaml_toys="./overrides/toysForProfiling.yaml"

safePrior="${newPrior//./_}"

base="${input_yaml##*/}"; base="${base%.*}"
override_dir="./macros/profiling/overrides/"
configs_json_dir="./macros/profiling/configs_json/"
configs_injector_dir="./macros/profiling/injector/"

# ensure dirs
cd ../..
mkdir -p "$override_dir" "$configs_json_dir"

override_output_yaml="$override_dir/${paramName}_${safePrior}.yaml"
json_out="$configs_json_dir/config_DUNE_With_${paramName}_${safePrior}.json"

injector_yaml="$configs_injector_dir/parameterInjector.yaml"
modif_injector="$configs_injector_dir/Inj_${paramName}_${safePrior}.yaml"

awk -v N="$occurrence" '
  BEGIN { count=0 }
  # Flip the Nth isFrozen line only:
  /^[[:space:]]*isFrozen:/ && ++count==N {
    match($0, /^[[:space:]]*/)                # grab indent
    printf "%sisFrozen: true\n", substr($0, RSTART, RLENGTH)
    next                                       # skip default print
  }
  { print }                                    # everything else unchanged
' "$input_yaml" > "$override_output_yaml"

awk -v N="$occurrence" -v V="$newPrior" '
  BEGIN { idx=0; fixNext=0 }

  # Case A: inline item:  - { name: "X", value: ... }
  /^[[:space:]]*-[[:space:]]*{[[:space:]]*name:[[:space:]]*"[^"]+"[[:space:]]*,[[:space:]]*value:/ {
    idx++
    if (idx==N) {
      # keep name as-is, just replace whole line with new value (preserving indent)
      line=$0
      match(line,/^[[:space:]]*/); lead=substr(line,RSTART,RLENGTH)
      # extract the name field to reprint it intact
      name=""
      if (match(line,/name:[[:space:]]*"([^"]+)"/,m)) name=m[1]
      printf "%s- { name: \"%s\", value: %s }\n", lead, name, V
      next
    }
    print; next
  }
  { print }

' "$injector_yaml" > "$modif_injector"

#produce json
gundamConfigUnfolder \
  -c config_DUNE.yaml \
  -of "$override_output_yaml" \
  -o "$json_out"

gundamFitter -a -c config_DUNE.yaml -of "$override_output_yaml" --inject-parameters "$modif_injector" --out-dir ./outputs/profiling/ -t 8
